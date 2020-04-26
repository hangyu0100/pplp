/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
 * This class provides functions to solve parametric linear programming.
 * _constraints stores the tableau:
 *        var1  var2  ...  varn  slack1  slack2  ... slackn  constant
 *  row1  coef  coef       coef    1       0           0         c   
 *  row2  coef  coef       coef    0       1           0         c
 *  rown  coef  coef       coef    0       0           1         c
 * _objective stores the tableau:
 *        var1  var2  ...  varn  slack1  slack2  ... slackn  constant
 *  para1 coef  coef       coef    0       0           0        0
 *  para2 coef  coef       coef    0       0           0        0
*******************************************************************************/

#ifndef _RAYTRACING_PLP
#define _RAYTRACING_PLP

// begin of debug info
#ifdef DEBUGINFO_PLP
#include <mutex>
#include <iostream>
#endif
// end of debug info

#include <vector>
#include <list>
#include <deque>
#include <chrono>
#include <mutex>

#include "pplp_polyhedron.h"
#include "pplp_glpkInterface.h"
#include "pplp_raytracing.h"

#ifdef _TBB
#define VERIMAG_POLYHEDRA_PLP_TBB 1
#define VERIMAG_POLYHEDRA_TBB_CONCURRENT_STRUCTURES 1

#include "tbb/parallel_do.h"

#else // NO TBB
#include <set>
#include <map>

#ifdef _OPENMP
#define VERIMAG_POLYHEDRA_PLP_OPENMP
#endif
#endif

#include <atomic>
#include <list>
#include <cstdint>
#include <queue>

#include <map>
#include <unordered_set>
#include <unordered_map>

namespace std {
  template<typename A, typename B> struct hash<std::pair<A, B>> {
    std::hash<A> hash_A;
    std::hash<B> hash_B;
  public:
    size_t operator() (const std::pair<A, B>& item) const {
      return hash_A(item.first) * 0xDEAD1133BEEFUL + hash_B(item.second);
    }
  };

  template<> struct hash<RMatrix> {
    std::size_t operator()(const RMatrix& m) const {
      std::size_t r = 0;
      for(int i=0; i<m.rows(); i++) {
	for(int j=0; j<m.cols(); j++) {
	  r = (r * 0xDEADBABEUL) + (m.at(i,j).num() % 1456789063UL).to<ulong>();
	  r = (r * 0xDEADCABEUL) + (m.at(i,j).den() % 1456789067UL).to<ulong>();
	}
      }
      return r;
    }
  };

  template<> struct hash< std::vector<int> > {
    std::size_t operator()(const std::vector<int>& vec) const {
      std::size_t r = 0;
      for(unsigned i = 0; i < vec.size(); ++ i) {
        r = (r * 0xDEADBABEUL) + vec[i] % 1456789063UL ;
      }
      return r;
    }
  };

  template<> struct less<RMatrix> {
    bool operator()( const RMatrix& a, const RMatrix& b) const {
      
      if (a.rows() < b.rows()) return true;
      if (a.rows() > b.rows()) return false;
      assert(a.rows() == b.rows());

      if (a.cols() < b.cols()) return true;
      if (a.cols() > b.cols()) return false;
      assert(a.cols() == b.cols());

      // at this point a and b have same number of lines and columns

      for(int i=0; i<a.rows(); i++) {
	for(int j=0; j<a.cols(); j++) {
	  auto ca = a.at(i, j), cb = b.at(i, j);
	  if (ca < cb) return true;
	  if (ca > cb) return false;
	}
      }
      // equality
      assert(a == b);
      return false;
    }
  };
}

// DM todo
#ifdef NO_ATOMIC64
typedef std::atomic<std::uint32_t> atomic_counter;
#else
typedef std::atomic<std::uint64_t> atomic_counter;
#endif

extern atomic_counter spin_count, work_microseconds;

namespace PPLP {

struct Region {
  Region(const Polyhedron& r_f, RMatrix&& r_r, Point&& p, const std::vector<int>& b,
      std::vector<int>&& ird) : f(std::move(r_f)), r(std::move(r_r)),
      point(std::move(p)), basis(std::move(b)), ird_cons(std::move(ird)) {
    #ifdef PLP_CHECKER_ON
    found_adj.resize(ird_cons.size(), false) ;
    #endif
      }
  Polyhedron f ;
  RMatrix r ;
  Point point ; 
  std::vector<int> basis ;
  std::vector<int> ird_cons ;
  std::vector<bool> found_adj ;
} ;

// this structure contains the optimal functions and their regions
struct Optimal {
  Optimal(Region&& reg, RMatrix&& opt)
      : region(std::move(reg)), optimal(std::move(opt)) {} ;
  Region region ;
  RMatrix optimal ;
} ;

#ifdef VERIMAG_POLYHEDRA_TBB_CONCURRENT_STRUCTURES
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_unordered_map.h"

typedef tbb::concurrent_vector<Optimal> optimal_container ;
typedef tbb::concurrent_unordered_map< std::vector<int>, int,
  std::hash< std::vector<int> > > base_map;
typedef tbb::concurrent_unordered_set<std::vector<int>,
  std::hash< std::vector<int> > > base_set;

#define LOCKFREE_SEEN_BASES
#define LOCKFREE_REGIONS 1
#define LOCKFREE_REGIONS_WAIT_FOR_CONDITION 1
#else

#ifdef VERIMAG_POLYHEDRA_PLP_OPENMP
#include "simple_concurrent_vector.h"
typedef simple_concurrent_vector<Optimal> optimal_container;

#define LOCKFREE_REGIONS 1
#define LOCKFREE_REGIONS_WAIT_FOR_CONDITION 1
#else
typedef std::vector<Optimal> optimal_container;

#endif

typedef std::map< std::vector<int>, int> base_map;
typedef std::set<std::vector<int>> base_set;
#endif

#if LOCKFREE_REGIONS
#if LOCKFREE_REGIONS_WAIT_FOR_CONDITION
#include <mutex>
#include <condition_variable>
#elif defined(__x86_64) || defined(__i686)
#include <immintrin.h>
#define LOCKFREE_REGIONS_PAUSE 1
#endif
#else
#include <mutex>
#endif

struct RegionIdx {
  RegionIdx() { on_boundary = -1 ; null_region = false ; }
  RegionIdx(int index, bool null = false, int bd = -1) 
    : on_boundary(bd), null_region(null), idx(index) {}
  int on_boundary ;
  bool null_region ;
  int idx ;
} ;

struct Objective {
  Objective(int rowNum, int colNum)
      : f(rowNum, colNum), r(rowNum, colNum) {}
  Objective(const RMatrix& obj) : r(obj) {}
  Matrix f ;
  RMatrix r ;
} ;

struct Task {
  Task(const RegionIdx& region, int frontier_idx, Vector&& p)
  : from_region(region), frontier(frontier_idx), point(std::move(p)) {}
  Task(const RegionIdx& region, int frontier_idx, const Vector& p)
  : from_region(region), frontier(frontier_idx), point(p) {}

  RegionIdx from_region ;
  int frontier ;
  Vector point ;
} ;

typedef std::vector<Task> Worklist_t;

#ifdef VERIMAG_POLYHEDRA_PLP_TBB
typedef tbb::parallel_do_feeder<Task> task_feeder;
#else
class task_feeder {
  Worklist_t& worklist;
  
public:
  task_feeder(Worklist_t& worklist0) : worklist(worklist0) { };

  void add(const Task& task) {
    worklist.push_back(task) ;
  }
  void add(Task&& task) {
    worklist.push_back(task) ;
  }
};
#endif

class Plp {
public:
  friend class TbbOperator ;
  friend class TbbParallel ;
    
  Plp(Polyhedron& poly, Worklist_t& worklist,
      optimal_container* optimals,
      int nb_initial_points) ;
  // Plp constructor for projection
  Plp(Polyhedron& poly, int projNum, const std::vector<int>& idx, Worklist_t& worklist,
      optimal_container* optimals, int nb_initial_points) ;
  // Plp constructor for convex hull
  Plp(Polyhedron& poly1, Polyhedron& poly2, Worklist_t& worklist,
      optimal_container* optimals,
      int nb_initial_points) ;
    
  bool IsNullRegion(const RegionIdx& idx) const {
    return idx.null_region ;
  }

  bool PointOnBoundary(const RegionIdx& idx) const {
    return idx.on_boundary != -1 ;
  }

  bool IsConstantRegion(const RegionIdx& idx) const {
    return _optimals->at(idx.idx).optimal.is_zero() ;
  }

  RMatrix& get_optimal(int idx) const {
    return _optimals->at(idx).optimal ; 
  }

  RMatrix& get_region_r(int idx) const {
    return _optimals->at(idx).region.r ; 
  }

  Polyhedron& get_region_f(int idx) const {
    return _optimals->at(idx).region.f ; 
  }

  void set_region_adj(int rIdx, int cIdx) {
    _optimals->at(rIdx).region.found_adj[cIdx] = true ;
  }
  bool get_region_adj(int rIdx, int cIdx) {
    return _optimals->at(rIdx).region.found_adj[cIdx] ;
  }

  std::vector<int>& get_region_basis(int idx) const {
    return _optimals->at(idx).region.basis ; 
  }

  std::vector<int>& get_region_ird_cons(int idx) const {
    return _optimals->at(idx).region.ird_cons ; 
  }

  RMatrix get_constraints() const {
    return _constraints ;
  }

  Matrix get_objective_f() const {
    return _objective.f ;
  }

  RMatrix get_objective_r() const {
    return _objective.r ;
  }

  /* Checks if the point is covered by one of the already seen regions.
   */
  RegionIdx CheckCovered(const Point& point);
  // int CheckCovered2(const Point& point) const;

  static RegionIdx GetFromRegion(const Task& currTask) {
    return currTask.from_region ;
  }
  static int GetFromFrontier(const Task& currTask) {
    return currTask.frontier ;
  }
  static Vector GetTaskPoint(const Task& currTask) {
    return currTask.point ;
  }

#ifndef LOCKFREE_SEEN_BASES
  std::mutex _seen_bases_lock;
#endif
  bool BasisAddExist(const std::vector<int>& idx) {
    std::vector<int> sortedIdx(idx) ;
    std::sort(sortedIdx.begin(), sortedIdx.end()) ;
#ifndef LOCKFREE_SEEN_BASES
    std::lock_guard<std::mutex> guard(_seen_bases_lock);
#endif
    return ! _seen_bases.insert(sortedIdx).second ;
  }

  std::vector<int> get_proj_idx() const {
    return _proj_idx ;
  }
  
  unsigned get_optimals_count() const { return _optimals_count; }
  
  int get_constraint_num() const {
    return _constraint_num ;
  }
  int get_variable_num() const {
    return _variable_num ;
  }
  std::vector<int> get_non_negative_idx() const {
    return _non_negative_idx ;
  }

  bool AreAdjacent(const RegionIdx& from,
      int frontier, const RegionIdx& currIdx) ;

  static bool ConsRedundantR(const RMatrix& matrix, int idx) ;
  static std::vector<int> MinimizePolyR(const RMatrix& matrix) ;
  
  bool IsFlatRegionR(const RMatrix& matrix) ;

private:
  base_set _seen_bases;
  RMatrix _constraints ; 
  Objective _objective ;
  int _constraint_num ;
  int _variable_num ;
  int _parameter_num ;
  // _poly store the constraints of plp
  Polyhedron _poly ;
  Polyhedron _glpk_poly ;
  optimal_container* _optimals ;
  std::atomic<int> _optimals_count;
#if LOCKFREE_REGIONS
#if LOCKFREE_REGIONS_WAIT_FOR_CONDITION
  std::mutex _optimals_mutex;
  std::condition_variable _optimals_cv;
#endif
#else
  std::mutex _optimals_mutex;
#endif

  std::vector<int> _proj_idx ;
  Vector _normalization_point ;
  // this object is used in the Constructor, and it will NOT 
  // effects the parallel functions
  std::vector<int> _non_negative_idx ;
  //std::vector<int> PretreatConsF() ; 
  int StoreResult(const RMatrix& result, const std::vector<int>& basicIdx,
      const std::vector<int>& nonBasicIdx, bool testPara = false) ;
  void AddTasks(const Task& currTask, int optimalIndex,
      task_feeder& feeder) ;

  void AddPoint(const Point& currPoint, const RegionIdx& currRegionIdx,
      const RegionIdx& fromRegionIdx, int fromConsIdx, 
      task_feeder& feeder, bool newRegion) const ;
  //bool AreAdjacent(const RegionIdx& from,
    //  int frontier, const RegionIdx& currIdx) ;
  bool CheckAdjacent(const RegionIdx& from, const RegionIdx& currIdx,
      const RMatrix& point) const ;
  RMatrix Reconstruct(const std::vector<int>& basicIdx,
      const RMatrix& cons, const RMatrix& obj) const ;
  RMatrix Reconstruct(const std::vector<int>& basicIdx) const ;
  RNumber GetRational(double num, int scale) const ; 
  double CutDouble(double num, int scale) const ; 
  bool CheckFirstPoint(const Vector& central_point,
      const Vector& first_point) const ;
  
  void add_initial_points(const Vector& realPoint,
			  const int NB_INITIAL_POINTS,
			  Worklist_t& worklist);
  int GetScale(double bound, int dimension) ;
  Vector GetSimplePoint(const Vector& point, int scale) ;
  static std::vector<int> GetNonDupIdxByCol(const RMatrix& result, 
      const std::vector<int>& idx) ;
  static std::vector<int> GetNonDupIdxByRow(const RMatrix& result, 
      const std::vector<int>& idx) ;
  RMatrix GetSubMatrix(const RMatrix& ori, const std::vector<int>& idx) ;
  //void ShiftPoint(const Task& currTask, task_feeder& feeder) ;
  bool ChooseBasis() ;
  int GetBdIdx(const RMatrix& matrix, int rowIdx) ;
  bool MinimizeLp() ;
  Vector GetSimplexObj(const Vector& currPoint) ;
  RMatrix ReconstructCons(const std::vector<int>& basicIdx) const ;
  RMatrix ReconstructCons(const std::vector<int>& basicIdx,
      const RMatrix& cons, bool useDefSize = false) const ;

  // TODO
  //void ProjSubEqs() ;

  bool VerifyGlpkResult(const std::vector<int>& basis) ;
#ifdef VERIFY_PLP
  bool VerifyCentralPoint(const RMatrix& point, const Polyhedron& poly) ;
#endif

} ;

}

#endif
