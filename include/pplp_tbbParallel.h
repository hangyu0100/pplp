/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
 //TODO description
*******************************************************************************/

#ifndef _RAYTRACING_TBB
#define _RAYTRACING_TBB

#ifdef _TBB
#include "tbb/parallel_do.h"
#endif

#include "pplp_plp.h"
#include "pplp_simplex.h"


#if defined(DEBUGINFO_PLP) || defined(PRINT_WARNING)
#include <mutex>
#include <iostream>
extern std::mutex log_mtx ;
#endif


typedef std::pair<int,int> FrontierIdx ;

namespace PPLP {
class TbbParallel {
public:
  TbbParallel(Polyhedron& poly) : _poly(poly) {}
  // for convex hull
  TbbParallel(Polyhedron& poly1, Polyhedron& poly2) : _poly(poly1), _poly2(poly2) {}
  void PlpParallel(int nb_initial_points) ; 
  void PlpParallel( int projNum, int nb_initial_points,
      const std::vector<int>& idx = std::vector<int>() ) ; 
  void PlpConvexHull(int nb_initial_points) ;
  void Operation(Plp& plpSolver, Worklist_t& worklist, optimal_container& optimals) ; 
  Polyhedron GetOptPoly() ;

  std::vector<int> get_proj_idx() {
    return _proj_idx ;
  }
  
  RMatrix GetOptimalMatrix() ;

private:
  std::unordered_set<RMatrix> _optimals ;
  Polyhedron _poly ;
  Polyhedron _poly2 ;
  std::vector<int> _proj_idx ;

  //std::queue<FrontierIdx> _to_find_adj ;
  void SearchMissedRegion(Plp& plp, std::queue< std::pair<int,int> >& check,
      const FrontierIdx& curr) ;
} ;

class TbbOperator {
public:
  TbbOperator(Plp* plpSolver, bool proj) : _plp(plpSolver), _projection(proj) {} ;
  void operator()(const Task& currTask, 
      task_feeder& feeder) const;
private:
  Plp* _plp ;
  bool _projection ;
} ;

}

#endif
