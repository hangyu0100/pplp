/*******************************************************************************
 * Copyright (c) 2018 Sep. Verimag. All rights reserved.
 * @author Hang YU
 * This class provides an algorithm of Simplex.
 * The constraints are in rational numbers, and the objective function is in
 * floating point numbers.
 * To solve a feasibility problem, the objective function is not needed. 
*******************************************************************************/

#ifndef _RAYTRACING_SIMPLEX
#define _RAYTRACING_SIMPLEX

#include <vector>
#include "pplp_polyhedron.h"
#include "pplp_raytracing.h"

#ifdef DEBUGINFO_SIMPLEX
#include <mutex>
#include <iostream>
#endif

namespace PPLP {

class Simplex {
public:
  Simplex( const RMatrix& cons, const RMatrix& obj, int consNum, int variNum,
      const Vector& point=Vector() ) ;
  Simplex( const RMatrix& cons, int consNum, int variNum,
      const Vector& point=Vector() ) ;
  std::vector<int> Solve(bool getFeasible=false) ;
  Vector GetFinalObjective() ;
  RMatrix get_cons_feasible() {
    return _cons_feasible ;
  }
  RMatrix get_feasible_vertex() {
    return _feasible_vertex ;
  }

private:  
  RMatrix PivotCons(const RMatrix& cons, int rowIdx, int colIdx) ;
  RMatrix Reconstruct(const std::vector<int>& basicIdx) const ;
  RMatrix Reconstruct(const std::vector<int>& basicIdx,
      const RMatrix& cons, const RMatrix& obj) const ;
  int GetOutBasicIdx(int inBasicCol, const RMatrix& coeff,
      const std::vector<int>& basicIdx) ;
  int IsOptimized(const RMatrix& objR, const std::vector<int>& nonBasicIdx) ;
  RMatrix _cons ;
  RMatrix _obj ;
  int _constraint_num ;
  int _variable_num ;
  // use Vector currently for testing
  // TODO: use template to adapt both Vector and RMatrix as obj 
  Vector _point ;
  RMatrix _cons_feasible ;
  RMatrix _feasible_vertex ;
} ;

}

#endif
