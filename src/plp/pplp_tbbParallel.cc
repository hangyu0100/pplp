#include "pplp_tbbParallel.h"

#include <set>

using namespace PPLP ;

#define PRINT_THREAD_ID 0
#if PRINT_THREAD_ID
#include <iostream>
#endif

// DM TODO
atomic_counter work_microseconds;
class AccumulateTime {
  typedef std::chrono::steady_clock clock;
  typedef std::chrono::time_point<clock> time_point;
  
  time_point start;

public:
  AccumulateTime() : start(clock::now()) { }
  ~AccumulateTime() {
    time_point end(clock::now());
    uint64_t delta =  std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    work_microseconds += delta;
  }
};

#if defined(VERIMAG_POLYHEDRA_PLP_TBB)
#elif defined (VERIMAG_POLYHEDRA_PLP_OPENMP)

void openmp_do_internal(const Worklist_t& worklist, const TbbOperator& operation);

void openmp_do_task(Task task, const TbbOperator& operation) {
  Worklist_t worklist;
  task_feeder feeder(worklist);
  operation(task, feeder);
  openmp_do_internal(worklist, operation);
}

void openmp_do_internal(const Worklist_t& worklist, const TbbOperator& operation) {
  for (Task task : worklist) {
#pragma omp task untied mergeable shared(operation)
    {
      openmp_do_task(std::move(task), operation);
    }
  }
}

void openmp_do(const Worklist_t& worklist, const TbbOperator& operation) {
  #pragma omp parallel
  {
    #pragma omp master
    {
      openmp_do_internal(worklist, operation);
    }
    // DON'T DO THAT flint_cleanup();
  }
}

#else
void serial_do(Worklist_t& worklist, const TbbOperator& operation) {
  task_feeder feeder(worklist);
  
  while(! worklist.empty()) {
    auto it = worklist.begin();
    Task task(std::move(*it));
    worklist.erase(it);
    operation(task, feeder);
  }
}
#endif

// for minimization
void TbbParallel::PlpParallel(int nb_initial_points) {
  optimal_container optimals ;
  Worklist_t worklist;
  Plp plpSolver(_poly, worklist, &optimals, nb_initial_points) ;
  // std::cout << "initial worklist = " << worklist.size() << std::endl;

#if defined(VERIMAG_POLYHEDRA_PLP_TBB)
  tbb::parallel_do(worklist, TbbOperator(&plpSolver, false) ) ;
#elif defined(VERIMAG_POLYHEDRA_PLP_OPENMP)
  openmp_do(worklist, TbbOperator(&plpSolver, false) ) ;
#else
  serial_do(worklist, TbbOperator(&plpSolver, false) ) ;
#endif

  assert(optimals.size() == plpSolver.get_optimals_count());
  std::unordered_set<RMatrix> minimized_optimals ;
  for (unsigned i = 0; i < plpSolver.get_optimals_count(); ++ i) {
    minimized_optimals.insert( plpSolver.get_optimal(i) ) ;
  }
  //std::cout << "regions: " << regions.size() << " optimals: " << optimals.size() << " minimized_optimals: " << minimized_optimals.size() << std::endl;
  _optimals = std::move(minimized_optimals) ;
}

// for projection
void TbbParallel::PlpParallel(int projNum, int nb_initial_points,
    const std::vector<int>& idx) {

#ifdef DEBUGINFO_PLP
  std::cout << "start parallel plp for projection" << std::endl ;
#endif

  optimal_container optimals ;
  Worklist_t worklist;
  std::vector<int> projIdx ;
  if ( idx.empty() ) {
    for (int i = 0; i < projNum; ++ i) {
     projIdx.push_back(i);
    }
  }
  else {
    projIdx = idx ;
  }

  Plp plpSolver(_poly, projNum, projIdx, worklist, &optimals, nb_initial_points) ;
 
  _proj_idx = plpSolver.get_proj_idx() ;
  Operation(plpSolver, worklist, optimals) ;
}

// for convex hull
void TbbParallel::PlpConvexHull(int nb_initial_points) {

#ifdef DEBUGINFO_PLP
  std::cout << "start parallel plp for projection" << std::endl ;
#endif

  optimal_container optimals ;
  Worklist_t worklist;
  Plp plpSolver(_poly, _poly2, worklist, &optimals, nb_initial_points) ;
  _proj_idx = plpSolver.get_proj_idx() ;
  Operation(plpSolver, worklist, optimals) ;
}

void TbbParallel::Operation(Plp& plpSolver, Worklist_t& worklist, optimal_container& optimals) {
  // if no solution
  int tmpRow = plpSolver.get_objective_f().rows() ;
  int tmpCol = plpSolver.get_objective_f().cols() ;
  // TODO change this
  // currently we set the objective function as 0 if the constraints
  // are infeasible or the whole space, i.e. bottom or top
  if ( plpSolver.get_objective_f().block(0,0, tmpRow, tmpCol-1).isZero() ) {

#ifdef PRINT_RESULT
  std::cout << "regions: 0 minimized_optimals: 0" << std::endl;
#endif

    return ;
  }

#if defined(VERIMAG_POLYHEDRA_PLP_TBB)
#ifdef DEBUGINFO_PLP
  std::cout << "start tbb operator" << std::endl ;
#endif

  tbb::parallel_do(worklist, TbbOperator(&plpSolver, true) ) ;
#elif defined(VERIMAG_POLYHEDRA_PLP_OPENMP)
  openmp_do(worklist, TbbOperator(&plpSolver, true) ) ;
#else
  serial_do(worklist, TbbOperator(&plpSolver, true) ) ;
#endif
  
// checker
#ifdef PLP_CHECKER_ON
  std::queue< std::pair<int,int> > toCheck ;
  for (unsigned i = 0; i < optimals.size(); ++ i) {
    int consSize = optimals[i].region.r.cols() ;
    for (int j = 0; j < consSize; ++ j) {
      int curr = plpSolver.get_region_adj(i, j) ;
      if ( ! curr ) {
        std::pair<int,int> curr(i, j) ;
        toCheck.push(curr) ; 
      }
    }
  }
  
  bool found ; 
  while ( ! toCheck.empty() ) {
    auto currCheck = toCheck.front() ;
    toCheck.pop() ;
    found = false ;
    int r1 = currCheck.first ;
    int c1 = currCheck.second ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "search adjacency of: " << "r: " << r1 << " c: " << c1 << std::endl ;
  log_mtx.unlock() ;
#endif

    RegionIdx rIdx1(r1) ;
    for (unsigned j = 0; j < plpSolver.get_optimals_count(); ++ j) {
      RegionIdx rIdx2(j) ;
      bool adj = plpSolver.AreAdjacent(rIdx1, c1, rIdx2) ;
      if (adj) {
        plpSolver.set_region_adj(r1, c1) ;
        found = true ;
        break ;
      }
    }
    if ( ! found ) {
      FrontierIdx currIdx(r1,c1) ;
      SearchMissedRegion(plpSolver, toCheck, currIdx) ;
    }
  }

#endif

  assert(optimals.size() == plpSolver.get_optimals_count());
  std::set<RMatrix> minimized_optimals ;
  for (unsigned i = 0; i < plpSolver.get_optimals_count(); ++ i) {
    if ( ! plpSolver.IsConstantRegion(i) ) {
      minimized_optimals.insert( plpSolver.get_optimal(i) ) ;
    }
  }
  std::unordered_set<RMatrix> minimized_optimals_hashed ;
  for (unsigned i = 0; i < plpSolver.get_optimals_count(); ++ i) {
    if ( ! plpSolver.IsConstantRegion(i) ) {
      minimized_optimals_hashed.insert( plpSolver.get_optimal(i) ) ;
    }
  }
  assert(minimized_optimals.size() ==  minimized_optimals_hashed.size());

#ifdef PRINT_RESULT
  std::cout << "regions: " << optimals.size() << " minimized_optimals: " << minimized_optimals.size() << std::endl;

  //for (auto& curr : minimized_optimals) {
  //  print(curr.transpose()) ;
  //}
#endif

  _optimals = std::move(minimized_optimals_hashed) ;


}

void TbbParallel::SearchMissedRegion(Plp& plp,
    std::queue< std::pair<int,int> >& check, const FrontierIdx& curr) {

#ifdef DEBUGINFO_PLP
  std::cout << "SearchMissedRegion() starts." << std::endl ;
#endif

  int rIdx, fIdx ;
  rIdx = curr.first ;
  fIdx = curr.second ;
  std::vector<int> basis = plp.get_region_basis(rIdx) ;
  RMatrix coeff( plp.ReconstructCons(basis) ) ; 
  RMatrix result( plp.Reconstruct(basis) );

#ifdef DEBUGINFO_PLP
    std::cout << "curr matrix: " << std::endl ;
    print(coeff) ;
    print(result) ; 
#endif

  int constIdx = coeff.cols()-1 ;
  int enteringIdx = plp.get_region_ird_cons(rIdx)[fIdx] ;
  RNumber zero ;
  RNumber ratio, currRatio ;

  int optimalIndex = RegionState::flat ;
  while (optimalIndex == RegionState::flat) {
    int leavingRow = -1 ;
    for (int i = 0; i < coeff.rows(); ++ i) {
      if ( coeff.at(i,enteringIdx) > zero ) {
        if ( coeff.at(i,constIdx).is_zero() ) {
          return ;
        }
        if (leavingRow == -1) {
          leavingRow = i ;
          ratio = coeff.at(i,constIdx) / coeff.at(i,enteringIdx) ;
        }
        else {
          currRatio = coeff.at(i,constIdx) / coeff.at(i,enteringIdx) ;
          if (currRatio < ratio) {
            leavingRow = i ;
            ratio = coeff.at(i,constIdx) / coeff.at(i,enteringIdx) ;
          }
          else if ( currRatio == ratio && basis[i] > basis[leavingRow] ) {
            leavingRow = i ;
            ratio = coeff.at(i,constIdx) / coeff.at(i,enteringIdx) ;
          }
        }
      }
    }
    if (leavingRow != -1) {
      basis[leavingRow] = enteringIdx ;
      std::vector<int> currNonBasic ;
      for (int i = 0; i < plp.get_variable_num(); ++ i) {
        if ( std::find(basis.begin(), basis.end(), i) == basis.end() ) {
          currNonBasic.push_back(i) ;
        }
      }

#ifdef DEBUGINFO_PLP
      std::cout << "new basis: " ;
      for (int idx : basis) {
        std::cout << idx << " " ;
      }
      std::cout << std::endl ;
#endif


      
      if ( plp.BasisAddExist(basis) ) {

#ifdef DEBUGINFO_PLP
      std::cout << "curr basis exists" << std::endl ; 
#endif

        return ;
      }
      
      RMatrix tmpCoeff( plp.ReconstructCons(basis) ) ; 
      RMatrix tmpResult(plp.Reconstruct(basis));

#ifdef DEBUGINFO_PLP
      std::cout << "new matrix: " << std::endl ;
      print(tmpCoeff) ;
      print(tmpResult) ; 
#endif

      optimalIndex = plp.StoreResult(tmpResult, basis, currNonBasic) ;

      if (optimalIndex == RegionState::flat) {
        // leaving row didn't change, find out flat cons
        for (int j = 0; j < tmpResult.cols()-1; ++ j) {
          if (j == enteringIdx) continue ;
          if ( Tool::IsZeroCol(tmpResult, j) ) continue ;
          if ( Tool::ColumnsEq(tmpResult, enteringIdx, j, true) ) {
            enteringIdx = j ;
            break ;
          }
        }
        Tool::CopyMatrix(tmpCoeff, coeff) ;
        Tool::CopyMatrix(tmpResult, result) ;
      }
    } 

    // add new frontier to the list
    if (optimalIndex != RegionState::flat) {
      int consNum = plp.get_region_f(optimalIndex).get_constraint_num() ;
      for (int i = 0; i < consNum; ++ i) {
        std::pair<int,int> curr(optimalIndex, i) ;
        check.push(curr) ;
      } 
    }
  }
    
#ifdef DEBUGINFO_PLP
    if (optimalIndex == RegionState::flat) {
      std::cout << "Cannot find new face" << std::endl ;
    }
#endif

}

/*******************************************************************************
 * Minimize polyhedron with plp  
 * @return the minimized polyhedron
 * TODO maybe remove it later
*******************************************************************************/
Polyhedron TbbParallel::GetOptPoly() {
  int consNum = _optimals.size() ;
  if (consNum == 0) {
    return Polyhedron(0, 0) ;
  }

  int variNum = _optimals.begin()->rows()-1 ;
  Polyhedron newPoly(consNum, variNum) ;
  int i = 0 ;
  for (auto it = _optimals.begin(); it != _optimals.end(); ++ it) {
    // Transfer into double matrix
    Vector currVec = Tool::GetDoubleMatrix( *it ).col(0).transpose() ; 
    newPoly.SetConstraint( i, -currVec.head(variNum) ) ;
    newPoly.SetConstant( i, currVec(variNum) ) ;
    ++ i ;
  }
  return newPoly ;
}

RMatrix TbbParallel::GetOptimalMatrix() {
  int variNum = _poly.get_variable_num() ;
  RMatrix optMatrix(_optimals.size(), variNum+1) ;
  if (_optimals.size() == 0) {
    optMatrix.set_zero() ;
  }
  else {
    auto it = _optimals.begin() ;
    for (int i = 0; it != _optimals.end(); ++ it, ++ i) {
      for (int j = 0; j < variNum; ++ j) {
        optMatrix.at(i, j) = it->at(j, 0) ;
      } 
      optMatrix.at(i, variNum) = it->at(variNum, 0) ;
    }
  }
  return optMatrix ;
}



static void log_thread(const char* phase
#ifdef __GNUC__
		       __attribute__((unused))
#endif		       
		       ) {
#if PRINT_THREAD_ID
  static auto start_time = std::chrono::high_resolution_clock::now();
  auto current_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_time).count();
  {
    std::lock_guard<std::mutex> guard(print_mutex);
    std::cout << phase
	      << " id=" << std::this_thread::get_id()
	      << ", time=" << current_time
	      << std::endl;
  }
#endif
}

class EnterExit {
public:
  EnterExit() {
    log_thread("BEGIN");
  }
  ~EnterExit() {
    log_thread("END");
  }
};


void TbbOperator::operator()(const Task& currTask, task_feeder& feeder) const {
  AccumulateTime accumulator; 

#if SERIALIZE
  static std::mutex mutex;
  std::lock_guard<std::mutex> guard(mutex);
#endif
  EnterExit enter_exit;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "curr task: " << currTask.point << std::endl ;
  log_mtx.unlock() ;
#endif

  bool adjacent = true ;
  Point currPoint( Plp::GetTaskPoint(currTask).head(_plp->_parameter_num) ) ; 
  RegionIdx from = Plp::GetFromRegion(currTask) ;
  int frontier = Plp::GetFromFrontier(currTask) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "from " ;
  if ( ! _plp->IsNullRegion(from) ) { 
    if ( _plp->IsConstantRegion(from) ) std::cout << " constant " ;
    else std::cout << " normal " ;
    std::cout << "region: " << from.idx << " frontier: " << frontier << std::endl ;
  } else { std::cout << "null" << std::endl ; }
  log_mtx.unlock() ;
#endif

  RegionIdx currIdx = _plp->CheckCovered(currPoint);
  if ( ! _plp->IsNullRegion(currIdx) ) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "current point is covered by" ;
  if ( _plp->IsConstantRegion(currIdx) ) std::cout << " constant " ; 
  else std::cout << " normal " ; 
  std::cout << "region " << currIdx.idx ;
  if (_plp->PointOnBoundary(currIdx)) std::cout << " on boundary" << std::endl ;
  else std::cout << " inside" << std::endl ;
  log_mtx.unlock() ;
#endif

    if ( ! _plp->IsNullRegion(from) ) {
      adjacent = _plp->AreAdjacent(from, frontier, currIdx) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  if (_plp->IsConstantRegion(from)) std::cout << "constant region " ;
  else std::cout << "normal region " ;
  std::cout << from.idx << " and " ;
  if (_plp->IsConstantRegion(currIdx)) std::cout << "constant region " ;
  else std::cout << "normal region " ;
  std::cout << currIdx.idx << " are adjacent? " << adjacent << std::endl ;
  log_mtx.unlock() ;
#endif

      if (! adjacent) {
        _plp->AddPoint(currPoint, currIdx, from, frontier, feeder, false) ;
      }
    }
  } 
  else {
    log_thread("NEW REGION");

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "new region" << std::endl ;
  log_mtx.unlock() ;
#endif

    // if the point does not belong to the known regions, 
    // get the new region
    
    bool glpkRes ;
    GlpkInterface glpkInter ;
    if (_projection) {
      Vector glpkObj = _plp->GetSimplexObj( Plp::GetTaskPoint(currTask) ) ;
      glpkRes = glpkInter.Simplex( _plp->_glpk_poly,
          glpkObj, true, false, true, _plp->get_non_negative_idx() ) ;
    }
    else {
      Vector lpObj = Plp::GetTaskPoint(currTask) *
          _plp->_objective.f.block(0, 0, _plp->_parameter_num+1, _plp->_variable_num) ;
      glpkRes = glpkInter.Simplex(_plp->_poly, lpObj, true) ;
    }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "glpk simplex OK" << std::endl ;
  log_mtx.unlock() ;
#endif

    // no optimal solution
    if (glpkRes == false) {
      return ;
    }
    std::vector<int> basicIdx( glpkInter.get_basic_idx() ); 
    std::vector<int> nonBasicIdx( glpkInter.get_non_basic_idx() ); 
    // the last item of the column is -optimal
    // check whether this region (given by a basis) has already been seen
    
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "basis: " << std::endl ;
  for (auto& curr : basicIdx) {
    std::cout << curr << " " ;
  }
  std::cout << std::endl ;
  glpkInter.PrintValue() ;
  log_mtx.unlock() ;
#endif

    bool verifyGlpk = _plp->VerifyGlpkResult(basicIdx) ;

#ifdef ALWAYS_RATIONAL
    verifyGlpk = false ;
#endif

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Verify result of GLPK: " << verifyGlpk << std::endl ; 
  log_mtx.unlock() ;
#endif

    if ( ! verifyGlpk) {

#ifdef PRINT_WARNING
  log_mtx.lock() ;
  std::cout << "Warning: VerifyGlpk(): infeasible solution from GLPK. Fix by hand." << std::endl ;
  log_mtx.unlock() ;
#endif
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Warning: VerifyGlpk(): infeasible solution from GLPK. Fix by hand." << std::endl ;
  log_mtx.unlock() ;
#endif

      Simplex exactSolver(_plp->get_constraints(), _plp->get_objective_r(),
          _plp->get_constraint_num(), _plp->get_variable_num(), currTask.point) ;
      basicIdx = exactSolver.Solve() ;


      // no feasible or optimal solution
      if (basicIdx.size() == 0) {
        return ;
      }
      nonBasicIdx.clear() ; 
      for (int i = 0, k = 0; i < _plp->_variable_num; ++ i) {
        if ( i != basicIdx[k] ) {
          nonBasicIdx.push_back(i) ;
        }
        else {
          ++ k ;
        }
      }
    }

#ifndef PLP_CHECKER_ON
    if ( _plp->BasisAddExist(basicIdx) ) {
      return ;
    }

    RMatrix result(_plp->Reconstruct(basicIdx));
    int optimalIndex = _plp->StoreResult(result, basicIdx, nonBasicIdx) ;
    if (optimalIndex >= 0) {
      _plp->AddTasks(currTask, optimalIndex, feeder) ;
    }

#else


    
/*******************************tmp code********************************
 * If we will find some regions here, we do not need to check adjacency,
 * as we cross one frontier a time, so it must be adjacent
 * Besides, we do not need to check feasibility, as it won't change if
 * the dictionary minus ax=0.
***********************************************************************/
    // TODO let VerifyGlpk use this as the parameter directly
    RMatrix coeff( _plp->ReconstructCons(basicIdx) ) ; 
    int constIdx = coeff.cols() - 1 ;
    std::vector<int> degRows, nonDegRows ;
    for (unsigned i = 0; i < coeff.rows(); ++ i) {
      if ( coeff.at(i,constIdx).is_zero() ) {
        degRows.push_back(i) ;
      }
      else {
        nonDegRows.push_back(i) ;
      }
    }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "degenerate row pushed: " ;
    for (int idx : degRows) {
      std::cout << idx << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif

    // if no degeneracy
    if ( degRows.empty() ) {
      if ( _plp->BasisAddExist(basicIdx) ) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "curr basis exists" << std::endl ; 
  log_mtx.unlock() ;
#endif

        return ;
      }

      RMatrix result(_plp->Reconstruct(basicIdx));

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "reconstructed result:" << std::endl ;
  print(result) ;
  std::cout << "GLPK objective: " ;
  double total, tmp ;
  for (int j = 0; j < result.cols()-1; ++ j) {
    total = 0 ;
    for (int i = 0; i < result.rows(); ++ i) {
      tmp = result.at(i,j).num().to<double>() / result.at(i,j).den().to<double>() ;
      total += currTask.point(i) * tmp ; 
    }
    std::cout << total << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

      int optimalIndex = _plp->StoreResult(result, basicIdx, nonBasicIdx) ;
      if (optimalIndex >= 0) {
        _plp->AddTasks(currTask, optimalIndex, feeder) ;
      }

      return ;
    }
    // if there is degeneracy, compute all the regions share the same face

    unsigned chooseNum = degRows.size() ;
    std::vector<int> chooseFirst ;
    for (int idx : degRows) {
      chooseFirst.push_back( basicIdx[idx] ) ;
    } 
    chooseFirst.insert( chooseFirst.end(), nonBasicIdx.begin(), nonBasicIdx.end() ) ;
    std::sort( chooseFirst.begin(), chooseFirst.end() ) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "choose first from: " ;
  for (int idx : chooseFirst) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

    std::vector< std::vector<int> > basisLexico ;
    Tool::BinomialCoeff(chooseFirst, chooseNum, basisLexico) ;
    
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "all the basis in lexico order: " ;
  for (const auto& vec : basisLexico) {
    for (int idx : vec) {
      std::cout << idx << " " ;
    }
    std::cout << std::endl ;
  }
  log_mtx.unlock() ;
#endif

    std::vector<int> firstBasis ;
    bool found = false ;
    for (unsigned i = 0; i < basisLexico.size(); ++ i) {
      firstBasis = basisLexico[i] ;
      for (int idx : nonDegRows) {
        firstBasis.push_back( basicIdx[idx] ) ;
      }

      
#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "first basis: " ;
    for (int idx : firstBasis) {
      std::cout << idx << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif

      RMatrix firstCoeff( _plp->ReconstructCons(firstBasis) ) ; 
      if ( firstCoeff.is_empty() ) {
        continue ;
      }

      std::vector<int> firstNonBasic ;
      for (int i = 0; i < _plp->get_variable_num(); ++ i) {
        if ( std::find(firstBasis.begin(), firstBasis.end(), i) == firstBasis.end() ) {
          firstNonBasic.push_back(i) ;
        }
      }

      RMatrix tmpResult(_plp->Reconstruct(firstBasis));
      // the result will not be stored
      // just check if the corresponding region is flat

      int state = _plp->StoreResult(tmpResult, firstBasis, firstNonBasic, true) ;

      // if found the first basis
      if ( state == RegionState::normal ) {
        found = true ;
        break ;
      }
    }

    if (! found) {
      return ;
    }

    std::queue< std::vector<ChooseFrom> > basicIdxList ;
    int consNum = _plp->get_constraint_num() ;
    RMatrix initPertMatrix(consNum, consNum) ;
    // make the variable with smallest subscript as the lexico smallest
    initPertMatrix.set_zero() ;
    for (int i = 0; i < consNum; ++ i) {
      initPertMatrix.at(i, consNum-i-1).set_one() ;
    }
    ChooseFrom chosenFirst( std::move(firstBasis), std::move(initPertMatrix), ConstraintIdx() ) ;
    std::vector< ChooseFrom > first ;
    first.push_back( std::move(chosenFirst) ) ;
    basicIdxList.push( std::move(first) ) ; 
    
    while( ! basicIdxList.empty() ) {
      std::vector< ChooseFrom > currList = basicIdxList.front() ;
      basicIdxList.pop() ;

      std::vector<int> entering ;
      std::map< int, std::vector<int> > dupIdx ;
      for (auto curr : currList) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "curr basis: " ; 
  for (auto idx : curr.currIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  
  std::cout << "pertMatrix: " ; print(curr.pertMatrix) ;
  log_mtx.unlock() ;
#endif 

        if ( _plp->BasisAddExist(curr.currIdx) ) {

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "curr basis exists" << std::endl ; 
    log_mtx.unlock() ;
#endif

          break ;
        }

        RMatrix coeff( _plp->ReconstructCons(curr.currIdx) ) ; 
        // if the current basis is not valid
        if ( coeff.is_empty() ) continue ;

        std::vector<int> degRows ;
        for (unsigned i = 0; i < coeff.rows(); ++ i) {
          if ( coeff.at(i,constIdx).is_zero() ) {
            degRows.push_back(i) ;
          }
        }

        std::vector<int> currNonBasic ;
        // find the non-basic variables
        for (int i = 0; i < _plp->get_variable_num(); ++ i) {
          //if ( std::find(curr.begin(), curr.end(), i) == curr.end() ) {
          if ( std::find(curr.currIdx.begin(), curr.currIdx.end(), i)
              == curr.currIdx.end() ) {
            currNonBasic.push_back(i) ;
          }
        }
        RMatrix tmpResult(_plp->Reconstruct(curr.currIdx));
        int regionState = _plp->StoreResult(tmpResult, curr.currIdx, currNonBasic) ;

#ifdef DEBUGINFO_PLP
        log_mtx.lock() ;
        std::cout << "constraints: " << std::endl ;
        print(coeff) ;
        std::cout << "result: " << std::endl ;
        print(tmpResult) ;
        std::cout << "degenerate row pushed: " ;
        for (int idx : degRows) {
          std::cout << idx << " " ;
        }
        std::cout << std::endl ;
        std::cout << "corresponding basis idx: " ;
        for (int idx : degRows) {
          std::cout << curr.currIdx[idx] << " " ;
        }
        std::cout << std::endl ;
        log_mtx.unlock() ;
#endif

        // begin to prepare the next pivot
        // if the region is flat
        bool flat = false ;
        bool foundFlat = false ;
        if ( regionState == RegionState::flat ) {
          flat = true ;
          int firstIdx, secondIdx ;
          for (unsigned i = 0; i < currNonBasic.size(); ++ i) {
            firstIdx = currNonBasic[i] ;
            for (unsigned j = i+1; j < currNonBasic.size(); ++ j) {
              secondIdx = currNonBasic[j] ;
              if ( Tool::IsZeroCol(tmpResult, firstIdx) ) {
                break ;
              }
              if ( Tool::IsZeroCol(tmpResult, secondIdx) ) {
                continue ;
              }
              if ( Tool::ColumnsEq(tmpResult, currNonBasic[i], currNonBasic[j], true) ) {
                entering.push_back( currNonBasic[i] ) ;
                entering.push_back( currNonBasic[j] ) ;
                foundFlat = true ;
                break ;
              }
            }
            if (foundFlat) {

#ifdef DEBUGINFO_PLP
      log_mtx.lock() ;
      std::cout << "Flat region found constraints: " ;
      for (int idx : entering) {
        std::cout << idx << " " ;
      }
      std::cout << std::endl ;
      log_mtx.unlock() ;
#endif
              break ;
            }
          }
        }
        // save information of adjacency
        
#ifdef PLP_CHECKER_ON
        else { 
          _plp->AddTasks(Task(RegionIdx(-1), -1, Vector()), regionState, feeder) ;
          entering = _plp->get_region_ird_cons(regionState) ;
          if (curr.from.regIdx != -1) {
            int regIdx = curr.from.regIdx ;
            int frontier = curr.from.frontierIdx ;
            for (unsigned i = 0; i < entering.size(); ++ i) {
              if( Tool::ColumnsEq(_plp->get_region_r(regIdx), tmpResult,
                  frontier, entering[i], true) ) {
                int lastOptIdx = _plp->get_optimals_count()-1 ;
                ConstraintIdx currIdx(lastOptIdx, i) ;
                _plp->set_region_adj(regIdx, frontier) ;
                _plp->set_region_adj(lastOptIdx, i) ;


#ifdef DEBUGINFO_PLP
      log_mtx.lock() ;
      std::cout << "found adjacency locally: " << regIdx << ": " << frontier
          << ", " << lastOptIdx << ": " << i << std::endl ;
      log_mtx.unlock() ;
#endif

                break ;
              }
            }
          }
        } 
#endif

#ifdef DEBUGINFO_PLP
      log_mtx.lock() ;
      if (flat && ! foundFlat) {
        std::cout << "Flat region cannot find constraints" << std::endl ;
        print(tmpResult.transpose()) ;
      }
      log_mtx.unlock() ;
#endif



        // for duplicated constraints
        for (unsigned i = 0; i < entering.size(); ++ i) {
          std::vector<int> currDupIdx ;
          for (unsigned j = 0; j < currNonBasic.size(); ++ j) {
            if ( entering[i] == currNonBasic[j] ) continue ;
            if ( Tool::IsZeroCol(tmpResult, currNonBasic[j]) ) continue ;
            if ( Tool::ColumnsEq( tmpResult, entering[i], currNonBasic[j] ) ) {
              currDupIdx.push_back( currNonBasic[j] ) ;
            }
          }
          dupIdx[ entering[i] ] = currDupIdx ;
        }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "entering: " << std::endl ;
    for (int idx : entering) {
      std::cout << idx << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif

        RNumber zero, currNum ;
        // lexico-smallest basis
        for (unsigned i = 0; i < entering.size(); ++ i) {
          auto it = dupIdx.find( entering[i] ) ;
          // all the columns duplicated
          std::vector<int> currIdxList ;
          if ( it != dupIdx.end() ) {
            currIdxList = it->second ;
          }
          currIdxList.insert( currIdxList.begin(), entering[i] ) ;
          std::vector<ChooseFrom> chooseFromList ;
          base_set choosed ;
          for (int currEntering : currIdxList) {
            int leavingIdx = -1 ;
            for (unsigned j = 0; j < degRows.size(); ++ j) {
              currNum = coeff.at(degRows[j],currEntering) ;
              if ( currNum > zero && Tool::LexicoPositive(curr.pertMatrix,degRows[j]) ) {
                if (leavingIdx == -1) {
                  leavingIdx = degRows[j] ;
                }
                else {
                  RMatrix vec1(1, consNum) ;
                  RMatrix vec2(1, consNum) ;
                  for (int k = 0; k < consNum; ++ k) {
                    vec1.at(0,k) = curr.pertMatrix.at(degRows[j],k) / currNum ;
                    vec2.at(0,k) = curr.pertMatrix.at(leavingIdx,k) / coeff.at(leavingIdx, currEntering) ;
                  }
                  if ( Tool::LexicoLess(vec1, vec2) ) {
                    leavingIdx = degRows[j] ;
                  }
                } 
              }
            }
            if (leavingIdx != -1) {
              std::vector<int> currBasicIdx(curr.currIdx) ;
              currBasicIdx[leavingIdx] = currEntering ;
              RMatrix newPertMatrix(consNum, consNum) ;
              RNumber den( coeff.at(leavingIdx, currEntering) ) ;
              for (int j = 0; j < consNum; ++ j) {
                newPertMatrix.at(leavingIdx,j) = curr.pertMatrix.at(leavingIdx,j) / den ;
              }
              for (int i = 0; i < consNum; ++ i) {
                if (i == leavingIdx) continue ;
                RNumber currCoeff( coeff.at(i, currEntering) ) ;
                for (int j = 0; j < consNum; ++ j) {
                  newPertMatrix.at(i,j) = curr.pertMatrix.at(i,j) -
                      currCoeff * newPertMatrix.at(leavingIdx,j) ;
                }
              }

              ConstraintIdx from(_plp->_optimals_count-1, i) ;
              if (flat) {
                from = curr.from ;
              }
              ChooseFrom newChoice(std::move(currBasicIdx), std::move(newPertMatrix), std::move(from)) ;
              chooseFromList.push_back(std::move(newChoice)) ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "chooseFrom pushed: " << std::endl ;
    for (int idx : currBasicIdx) {
      std::cout << idx << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif

            }

          }
          std::sort(chooseFromList.begin(), chooseFromList.end(), Tool::CmpFromList) ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "chooseFrom: " << std::endl ;
    for (const auto& currC : chooseFromList) {
      for (int c : currC.currIdx) {
        std::cout << c << " " ;
      }
      std::cout << std::endl ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif



          basicIdxList.push( std::move(chooseFromList) ) ;
        }

        break ;
      }
    }

#endif

/*******************************tmp code ends***************************/

  }
}

