#include "pplp_simplex.h"

using namespace PPLP ;

#ifdef DEBUGINFO_SIMPLEX
std::mutex log_mtx_spx;
#endif

Simplex::Simplex(const RMatrix& cons, const RMatrix& obj, int consNum, int variNum,
    const Vector& point) 
  : _cons(cons), _obj(obj), _constraint_num(consNum), _variable_num(variNum),
    _point(point), _cons_feasible(_constraint_num, _variable_num+1),
    _feasible_vertex(1, consNum) {
      _cons_feasible.set_zero() ;
    }

Simplex::Simplex(const RMatrix& cons, int consNum, int variNum, const Vector& point)
  : _cons(cons), _obj( RMatrix(0,0) ), _constraint_num(consNum),
    _variable_num(variNum), _point(point),
    _cons_feasible(_constraint_num, _variable_num+1), _feasible_vertex(1, consNum) {
      _cons_feasible.set_zero() ;
    }

std::vector<int> Simplex::Solve(bool getFeasible) {
  int nonBasicNum = _variable_num - _constraint_num ; 
  int constIdx = _cons.cols() - 1 ;
  std::vector<int> basicIdx ;
  std::vector<int> nonBasicIdx ;
  // get initial basis
  for (int i = 0; i < nonBasicNum; ++ i) {
    nonBasicIdx.push_back(i) ;
  }
  for (int i = nonBasicNum; i < _variable_num; ++ i) {
    basicIdx.push_back(i) ;
  }
  // check feasibility of initial directory
  bool oriFeasible = true ;
  RNumber zero ;
  for (int i = 0; i < _constraint_num; ++ i) {
    if (_cons.at(i, constIdx) < zero) {
      oriFeasible = false ;
      break ; 
    }
  } 
  
#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): check initial feasibility: " << oriFeasible << std::endl ; 
  std::cout << "Initial basic idx: " ;
  for (auto idx : basicIdx) {
    std::cout << idx << " " ;
  }
  std::cout << " Initial nonbasic idx: " ;
  for (auto idx : nonBasicIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  log_mtx_spx.unlock() ;
#endif

  // feasiblility phase
  bool optimized ;
  RMatrix cons_feasible_aux(_constraint_num, _variable_num+2) ;
  if ( ! oriFeasible) {
    RMatrix obj_auxiliary(1, _variable_num+2) ;
    obj_auxiliary.set_zero(); 
    obj_auxiliary.at(0, _variable_num) = 1 ;
    RMatrix cons_auxiliary(_constraint_num, _variable_num+2) ;
    int outIdx = -1 ;
    RNumber most, curr ;
    for (int i = 0; i < _constraint_num; ++ i) {
      for (int j = 0; j < _variable_num; ++ j) {
        cons_auxiliary.at(i, j) = _cons.at(i, j) ;
      }
      cons_auxiliary.at(i, _variable_num) = -1 ;
      curr = _cons.at(i, constIdx) ;
      cons_auxiliary.at(i, _variable_num+1) = curr ;
      if (curr < zero && curr < most) {
        most = curr ;
        outIdx = i ;
      }
    }
    
    if (outIdx == -1) {
      std::cerr << "Error, result of checking feasibility is wrong." << std::endl ;
      std::terminate() ;
    }
    nonBasicIdx.push_back(_variable_num) ;

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "start to find feasibility" << std::endl ; 
  print(cons_auxiliary) ;
  print(obj_auxiliary) ;
  std::cout << "the variable into basis: " << nonBasicIdx[nonBasicNum] << std::endl ;
  std::cout << "the variable out of basis: " << basicIdx[outIdx] << std::endl ;
  log_mtx_spx.unlock() ;
#endif

    std::swap(basicIdx[outIdx], nonBasicIdx[nonBasicNum]) ;
    std::sort( nonBasicIdx.begin(), nonBasicIdx.end() ) ;
    std::vector<int> sortedIdx = basicIdx ;
    std::sort( sortedIdx.begin(), sortedIdx.end() ) ;
    RNumber optVal ;
    RMatrix coeff( PivotCons( cons_auxiliary, outIdx, basicIdx[outIdx] ) ) ;
    RMatrix objR( Reconstruct(sortedIdx, cons_auxiliary, obj_auxiliary) ) ;

    do {
      optimized = true ;
      // do not need to consider the constant
      // check optimization, all the coefficients are non-negtive
      int inBasic, outBasic = -1 ;
      for (int j = 0; j < nonBasicNum+1; ++ j) {
        if (objR.at( 0, nonBasicIdx[j] ) < zero) {
          optimized = false ;
          inBasic = j ;
          break ;
        }
      }  
      if (optimized) {

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "found feasible:" << std::endl ;
  std::cout << "coeff: " ; print(coeff) ;
  std::cout << "objR: " ; print(objR) ;
  log_mtx_spx.unlock() ;
#endif

        optVal = objR.at(0, _variable_num+1) ;
        cons_feasible_aux = coeff ;
        break ;
      } 
      
      // find the first smallest positive ratio
      int inBasicCol = nonBasicIdx[inBasic] ;
      outBasic = GetOutBasicIdx(inBasicCol, coeff, basicIdx) ; 
      //unbounded
      if (outBasic == -1) {

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): first phase unbounded" << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

        return std::vector<int>() ;
      }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "coeff: " ; print(coeff) ;
  std::cout << "obj: " ; print(objR) ;
  std::cout << "the variable into basis: " << nonBasicIdx[inBasic] << std::endl ;
  std::cout << "the variable out of basis: " << basicIdx[outBasic] << std::endl ;
  log_mtx_spx.unlock() ;
#endif

      // change the basis
      std::swap(basicIdx[outBasic], nonBasicIdx[inBasic]) ;
      // sort in ascending order
      std::sort( nonBasicIdx.begin(), nonBasicIdx.end() ) ;
      std::vector<int> sortedIdx = basicIdx ;
      std::sort( sortedIdx.begin(), sortedIdx.end() ) ; 
      coeff = PivotCons( coeff, outBasic, basicIdx[outBasic] ) ;
      objR = Reconstruct(sortedIdx, coeff, objR) ;
    } while ( ! optimized) ;

// infeasible if the optimized value is not 0
    if ( ! optVal.is_zero() ) {

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): no feasible solution" << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

      return std::vector<int>() ;
    }
    else {
      // if it is degenerate
      if ( objR.at(0, _variable_num).is_zero() ) {

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): feasible solution degenerate" << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

        int outtmp = -1, intmp = -1 ;
        for (unsigned k = 0; k < basicIdx.size(); ++ k) {
          if (basicIdx[k] == _variable_num) {
            outtmp = k ;
            break ;
          }
        }
        for (unsigned k = 0; k < nonBasicIdx.size(); ++ k) {
          if ( ! objR.at( 0, nonBasicIdx[k] ).is_zero() ) {
            intmp = k ; 
            break ;
          }
        }
        if (outtmp == -1 || intmp == -1) {
          std::cerr << "Simplex(): first phase degeneracy checking error" << std::endl ;
          std::terminate() ;
        }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "variable into basis: " << nonBasicIdx[intmp] << std::endl ; 
  std::cout << "variable out of basis: " << basicIdx[outtmp] << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

        std::swap( basicIdx[outtmp], nonBasicIdx[intmp] ) ;
        cons_feasible_aux = PivotCons( coeff, outtmp, basicIdx[outtmp] ) ;
        std::sort( nonBasicIdx.begin(), nonBasicIdx.end() ) ;

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): feasible solution degenerate fixed" << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

      }
    }
    nonBasicIdx.erase(nonBasicIdx.end()-1) ;
  }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): found feasible: " << optimized << std::endl ; 
  if (optimized) {
    std::cout << "basic variable index: " ;
    for (unsigned i = 0; i < basicIdx.size(); ++ i) {
      std::cout << basicIdx[i] << " " ;
    }
    std::cout << std::endl ;
    std::cout << "nonbasic variable index: " ;
    for (unsigned i = 0; i < nonBasicIdx.size(); ++ i) {
      std::cout << nonBasicIdx[i] << " " ;
    }
    std::cout << std::endl ;
    if ( ! oriFeasible ) {
      std::cout << "auxiliary feasible constraints: " ; print(cons_feasible_aux) ;
    }
  }
  log_mtx_spx.unlock() ;
#endif

  if (oriFeasible) {
    for (int i = 0; i < _constraint_num; ++ i) {
      for (int j = 0; j < _variable_num; ++ j) {
        _cons_feasible.at(i, j) = _cons.at(i, j) ;
      }
      _cons_feasible.at(i, _variable_num) = _cons.at(i, constIdx) ;
    }
  }
  else {
    for (int i = 0; i < _constraint_num; ++ i) {
      for (int j = 0; j < _variable_num; ++ j) {
        _cons_feasible.at(i, j) = cons_feasible_aux.at(i, j) ; 
      }
      _cons_feasible.at(i, _variable_num) = cons_feasible_aux.at(i, _variable_num+1) ;
    }
  }
  for (unsigned i = 0; i < basicIdx.size(); ++ i) {
    _feasible_vertex.at(0,i) = _cons_feasible.at(i, _variable_num) ;
  } 

  if (getFeasible) {
    return basicIdx ;
  }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "start optimization phase" << std::endl ;
  log_mtx_spx.unlock() ;
#endif

  // start optimization phase
  RMatrix coeff(_cons_feasible) ;
  std::vector<int> sortedIdx = basicIdx ;
  std::sort( sortedIdx.begin(), sortedIdx.end() ) ;
  RMatrix objR( Reconstruct(sortedIdx) ) ;
  optimized = false ;
  do {
    int outBasic = -1 ;
    int inBasic = IsOptimized(objR, nonBasicIdx) ; 
    if (inBasic == -1) {
      optimized = true ;
    }
    if (optimized) {
      break ;
    }

    // find the first smallest positive ratio
    int inBasicCol = nonBasicIdx[inBasic] ;
    outBasic = GetOutBasicIdx(inBasicCol, coeff, basicIdx) ;
    //unbounded
    if (outBasic == -1) {

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Plp Simplex(): second phase unbounded" << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

      return std::vector<int>() ;
    }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "coeff: " ; print(coeff) ;
  std::cout << "obj: " ; print(objR) ;
  std::cout << "the variable into basis: " << nonBasicIdx[inBasic] << std::endl ;
  std::cout << "the variable out of basis: " << basicIdx[outBasic] << std::endl ;
  log_mtx_spx.unlock() ;
#endif
  
    // change the basis
    std::swap(basicIdx[outBasic], nonBasicIdx[inBasic]) ;
    // sort in ascending order
    std::sort( nonBasicIdx.begin(), nonBasicIdx.end() ) ;
    std::vector<int> sortedIdx = basicIdx ;
    std::sort( sortedIdx.begin(), sortedIdx.end() ) ;
    coeff = PivotCons( coeff, outBasic, basicIdx[outBasic] ) ;
    objR = Reconstruct(sortedIdx, coeff, objR) ;
  } while ( ! optimized) ;
  
  return basicIdx ;
}

RMatrix Simplex::PivotCons(const RMatrix& cons, int rowIdx, int colIdx) {
  int rowNum = cons.rows() ;
  int colNum = cons.cols() ;
  RMatrix newCons(cons) ;
  for (int j = 0; j < colNum; ++ j) {
    newCons.at(rowIdx, j) = cons.at(rowIdx, j) / cons.at(rowIdx, colIdx) ; 
  }
  RNumber ratio ;
  for (int i = 0; i < rowNum; ++ i) {
    if (i == rowIdx) continue ;
    ratio = newCons.at(i, colIdx) ; 
    for (int j = 0; j < colNum; ++ j) {
      newCons.at(i, j) = newCons.at(i, j) - ratio * newCons.at(rowIdx, j) ; 
    }
  }
  return newCons ;
}

/*******************************************************************************
 * Reconstructs the simplex tableau
 * @para constratins the tableau of constraints
 * @para obj the tableau of objective
 * @return the reconstructed tableau
*******************************************************************************/
RMatrix Simplex::Reconstruct(const std::vector<int>& basicIdx) const {
  // avoid zero rows in cons
  if (_cons.rows() != _constraint_num) {
    RMatrix newCons( Tool::GetBlock( _cons, 0, 0, _constraint_num, _cons.cols() ) ) ;
    return Reconstruct(basicIdx, newCons, _obj) ;
  }
  else {
    return Reconstruct(basicIdx, _cons, _obj) ;
  }
}

RMatrix Simplex::Reconstruct(const std::vector<int>& basicIdx, const RMatrix& cons, 
      const RMatrix& obj) const {
  int basicNum = basicIdx.size() ;
  RMatrix objMatrix(obj) ;
  int consNum = cons.rows() ;
  int variNum = cons.cols() - 1 ;
  RMatrix coeff(consNum, basicNum) ;
  coeff.set_zero() ;
  RMatrix objCoeff(objMatrix.rows(), basicNum) ;
  objCoeff.set_zero() ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < basicNum; ++ j) {
      coeff.at(i,j) = cons.at(i, basicIdx[j]) ;
    }
  }
  for (int i = 0; i < objMatrix.rows(); ++ i) {
    for (int j = 0; j < basicNum; ++ j) {
      // if glpk chose an row variable as a basic variable, the column in
      // objective should be 0
      if (basicIdx[j] <= variNum) {
        objCoeff.at(i,j) = objMatrix.at(i, basicIdx[j]) ;
      }
    }
  }

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "solve linear equation: " ;
  std::cout << " Lambda * " ; print(coeff) ;
  std::cout << " = " ; print(objCoeff) ;
  log_mtx_spx.unlock() ;
#endif

  // dixon vs fraction_free
  //fraction-free seems to have better performance
  
  auto lambda( coeff.transpose().solve_fraction_free( objCoeff.transpose() ) ) ;
  // avoid zero rows in cons
  /*
  if (lambda.rows() != cons.rows()) {
    RMatrix newCons( Tool::GetBlock( cons, 0, 0, lambda.rows(), cons.cols() ) ) ;
    return RMatrix(objMatrix - lambda.transpose() * newCons) ; 
  }
  else {
    return RMatrix(objMatrix - lambda.transpose() * cons) ; 
  } 
  */
  return RMatrix(objMatrix - lambda.transpose() * cons) ; 
}

int Simplex::GetOutBasicIdx(int inBasicCol, const RMatrix& coeff,
    const std::vector<int>& basicIdx) {
  int colNum = coeff.cols() ;
  RNumber ratio, curr, zero, one ;
  one.set_one() ;
  ratio = - one ;
  int outBasic = -1 ;
  // find the first smallest positive ratio
  for (int i = 0; i < _constraint_num; ++ i) {
    if (coeff.at(i, inBasicCol) <= zero) continue ;
    curr = coeff.at(i, colNum-1) / coeff.at(i, inBasicCol) ;
    if (ratio == -one) {
      ratio = curr ;
      outBasic = i ;
      continue ;
    }
    if (curr < ratio) {
      ratio = curr ; 
      outBasic = i ;
    }
    // find the one with smallest subscript
    else if ( curr == ratio && basicIdx[i] < basicIdx[outBasic] ) {
      outBasic = i ; 
    }
  }
  return outBasic ;
}

/*******************************************************************************
 * Check if the simplex reaches optimization
 * @return -1 if the optimal solution is found, otherwise reutrns the index of
 * the non-basic variable which will enter to the basis
*******************************************************************************/
int Simplex::IsOptimized(const RMatrix& objR, const std::vector<int>& nonBasicIdx) {
  int nonBasicNum = nonBasicIdx.size() ;
  // do not need to consider the constant
  Matrix objTmp(_obj.rows(), nonBasicNum) ;
  objTmp.setZero() ;
  double num, den ;
  // do not need to consider the constant
  for (int i = 0; i < _obj.rows(); ++ i) {
    for (int j = 0; j < nonBasicNum; ++ j) {
      num = objR.at( i, nonBasicIdx[j] ).num().to<double>() ;
      den = objR.at( i, nonBasicIdx[j] ).den().to<double>() ;
      objTmp(i, j) = num/den ;
    }
  } 
  Vector objF = _point * objTmp ;

#ifdef DEBUGINFO_SIMPLEX
  log_mtx_spx.lock() ;
  std::cout << "Simplex(): curr obj: " << objF << std::endl ; 
  log_mtx_spx.unlock() ;
#endif

  // check optimization, all the coefficients are non-negtive
  int inBasic = -1 ;
  double epsilon ;
  for (int j = 0; j < nonBasicNum; ++ j) {
    epsilon = std::numeric_limits<double>().epsilon() ;
    if ( Double::IsLessThan(objF(j) , 0, epsilon) ) {
      inBasic = j ;
      break ;
    } 
  }
  return inBasic ;
}
