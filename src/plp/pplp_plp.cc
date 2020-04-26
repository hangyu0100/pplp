/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
*******************************************************************************/
#include <utility>
#include <fstream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <random>
#include <iomanip>

#include "pplp_plp.h"
#include "pplp_double.h"
#include "pplp_apronInterface.h"
#include "pplp_simplex.h"

using namespace PPLP ;

#ifdef DEBUGINFO_PLP
std::mutex log_mtx;
#endif


// DM
atomic_counter spin_count;

/* DM: TODO
maybe make it deterministic?
*/

void Plp::add_initial_points(const Vector& realPoint,
			     const int NB_INITIAL_POINTS,
			     Worklist_t& worklist){

  // add the first points to the point stack
  Vector firstPoint(_parameter_num+1) ;
  RegionIdx nullRegion(-1, true) ;
  if (NB_INITIAL_POINTS == 0) {
    firstPoint.head(_parameter_num).setZero();
    firstPoint(_parameter_num) = 1.0;
    //worklist.emplace( FromIndex(-1,0), std::move(firstPoint) ) ;
    worklist.emplace_back(nullRegion, 0, std::move(firstPoint) ) ;
  } else {
    for(int k=0; k<NB_INITIAL_POINTS; k++) {
      do {
	firstPoint.head(_parameter_num).setRandom() ;
	firstPoint(_parameter_num) = 1.0;
      } while (firstPoint.head(_parameter_num) == realPoint)  ;
      //worklist.emplace( FromIndex(-1,k), std::move(firstPoint) ) ;
      worklist.emplace_back(nullRegion, k, std::move(firstPoint) ) ;
    }
  }
#if 0
  for(auto point : worklist) {
    std::cout << "direction: " << point.second.transpose() << std::endl;
  }
  std::cout << "end directions" << std::endl;
#endif
}

/*******************************************************************************
 * Construct plp for minimization
 * @para poly the polyhedron to be minimized
*******************************************************************************/
Plp::Plp(Polyhedron& poly,
	 Worklist_t& worklist,
         optimal_container* optimals,
         const int NB_INITIAL_POINTS)
    : _constraints(1, poly.get_constraint_num()+2),
      _objective(poly.get_variable_num()+1, poly.get_constraint_num()+2),
      _poly(1, poly.get_constraint_num()+1), _optimals(optimals),
      _optimals_count(0) {
  _constraint_num = 1 ;
  // lambda_i C_i + lambda_0 (for constant)
  _variable_num = poly.get_constraint_num() + 1;
  _parameter_num = poly.get_variable_num() ;
  // column num is _variable_num + 1, 
  // _parameter_num+1 as the last row is the constant
  int colNum = _variable_num + 1 ;
  _objective.f.setZero() ; 
  // construct the objective function according to the input polyhedron
  _objective.f.block(0, 0, _parameter_num, _variable_num-1) 
      = -poly.get_coefficients().transpose() ; 
  // Cx+b >= 0
  _objective.f.row(_parameter_num).head(_variable_num-1) = poly.get_constants() ;
  _objective.f(_parameter_num, _variable_num-1) = 1.0 ;

  // set the objective matrix of flint
  long int val ;
  _objective.r.set_zero() ;
  for (int i = 0; i < _parameter_num+1; ++ i) {
    for(int j = 0; j < colNum; ++ j) { 
      val = Tool::GetIntOfFloat( _objective.f(i, j) ) ;
      _objective.r.at(i, j) = val ;
    }
  }

  // construct the constraints
  _constraints.set_zero() ;
  // to normalization, evaluate the objective a point 
  // which is inside the input polyhedron
  Point tmpPoint = GlpkInterface::GetCentralPoint(poly) ;
  if ( tmpPoint.IsEmpty() ) {
    tmpPoint = GlpkInterface::GetSatPoint(poly) ;
  }
  if( tmpPoint.IsEmpty() ) {
    std::cerr << "Minimization failed. Cannot find a point inside the polyhedra" << std::endl ;
    std::terminate() ;
  }

  int cutScale = GetScale( poly.get_dis_lower_bound(), _parameter_num ) ;
  _normalization_point = GetSimplePoint( tmpPoint.get_coordinates(), cutScale) ;
  Vector equalCons(_variable_num) ;
  equalCons.head(_variable_num-1) =
      _normalization_point * -poly.get_coefficients().transpose()
      + poly.get_constants() ;
  equalCons(_variable_num-1) = 1.0 ;
  _poly.SetConstraint( 0, std::move(equalCons) ) ;
  _poly.SetConstant(0, 1.0) ;
  //_poly.SetOperator(0, ConstraintOperator::equal) ;

  // all constraints in _poly are Cx <= b
  // The matrix in glpk is:
  //       col1  col2  ...  coln - row1    row2  ...  rown <= constant
  // row1  coef  coef       coef    -1       0           0         c   
  // row2  coef  coef       coef    0       -1           0         c
  // rown  coef  coef       coef    0       0           -1         c
  // But we just consider the case that all the constraints are equalities,
  // so we do not add the row variables to our table _constraints
  RNumber curr ;
  for (int j = 0; j < _variable_num; ++ j) {
    curr = GetRational( equalCons(j), cutScale) ;
    _constraints.at(0, j) = std::move(curr) ;
  }
  // the constraint is Cx = 1
  _constraints.at(0, colNum-1) = 1 ;
  
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints: " << std::endl ; print(_constraints) ;
  std::cout << "objective: " << std::endl ; print(_objective.r) ;
  log_mtx.unlock() ;
#endif

  // add the first points to the point stack
  Vector firstPoint(_parameter_num+1) ;
  RegionIdx nullRegion(-1, true) ;
  if (NB_INITIAL_POINTS == 0) {
    firstPoint.head(_parameter_num).setZero();
    firstPoint(_parameter_num) = 1.0;
    //worklist.emplace( FromIndex(-1,0), std::move(firstPoint) ) ;
    worklist.emplace_back(nullRegion, 0, std::move(firstPoint) ) ;
  } else {
    for(int k=0; k<NB_INITIAL_POINTS; k++) {
      do {
	firstPoint.head(_parameter_num).setRandom() ;
	firstPoint(_parameter_num) = 1.0 ;
      } while (firstPoint.head(_parameter_num) == _normalization_point)  ;
      worklist.emplace_back(nullRegion, k, firstPoint) ;
    }
  }
}

/*******************************************************************************
 * Construct plp for projection
 * Randomly choose projNum variables to project
 * @para poly the polyhedron to be projected
*******************************************************************************/
Plp::Plp(Polyhedron& poly,
         int projNum,
         const std::vector<int>& idx,
	 Worklist_t& worklist,
         optimal_container* optimals,
         const int NB_INITIAL_POINTS)
    : _constraints(1+projNum, poly.get_constraint_num()+2),
      _objective(poly.get_variable_num()+1, poly.get_constraint_num()+2),
      _optimals(optimals), _optimals_count(0), _proj_idx(idx) {
  // lambda_i C_i + lambda_0 (for constant)
  //_variable_num = poly.get_constraint_num() + 1;
  _parameter_num = poly.get_variable_num() ;
  // column num is _variable_num + 1, 
  // _parameter_num+1 as the last row is the constant
  int initVariNum = poly.get_constraint_num() + 1 ;
  int colNum = initVariNum + 1 ;

  // TODO substitute eqs
  /*
  if (poly.get_eq_constraint_num() != 0) {
    ProjSubEqs() ;
  }
  */
 
  // construct the constraints
  _constraints.set_zero() ;
  // to normalization, evaluate the objective a point 
  // which is inside the input polyhedron

  Point tmpPoint = GlpkInterface::GetCentralPoint(poly) ;
  if ( tmpPoint.IsEmpty() ) {
    tmpPoint = GlpkInterface::GetSatPoint(poly) ;
  }
  if( tmpPoint.IsEmpty() ) {
    std::cerr << "Minimization failed. Cannot find a point inside the polyhedra" << std::endl ;
    std::terminate() ;
  }

  int cutScale = GetScale( poly.get_dis_lower_bound(), _parameter_num ) ;
  _normalization_point = GetSimplePoint( tmpPoint.get_coordinates(), cutScale) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Lower bound of diatance: " << poly.get_dis_lower_bound() << std::endl ;
  log_mtx.unlock() ;
#endif

#ifdef VERIFY_PLP
  RMatrix rationalPoint(_parameter_num, 1) ;
  for (int i = 0; i < _parameter_num; ++ i) {
    rationalPoint.at(i, 0) = GetRational(_normalization_point(i), cutScale) ;
  }
  bool pointVerified = VerifyCentralPoint(rationalPoint, poly) ;
  assert(pointVerified) ;
#endif

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "tmp central point OK: " << tmpPoint.get_coordinates().transpose()
      << std::endl ;
  std::cout << "real central point OK: " << _normalization_point.transpose()
      << std::endl ;
  log_mtx.unlock() ;
#endif

  Vector equalCons(initVariNum) ;
  equalCons.head(initVariNum-1) =
      _normalization_point * -poly.get_coefficients().transpose()
      + poly.get_constants() ;

  equalCons(initVariNum-1) = 1.0 ;
  std::vector<int> colIdx ;
  for (int i = 0; i < projNum; ++ i) {
    // avoid the constraints like 0lambda = 0
    if ( ! poly.get_coefficients().col(_proj_idx[i]).isZero() ) {
      colIdx.push_back( _proj_idx[i] ) ;
    }
  }
  int colSize = colIdx.size() ;
  // +1 for the normalization constraint
  int init_cons_num = colSize + 1 ;
  // all constraints are equalities
  _poly.SetSize(0, initVariNum, init_cons_num) ;
  _poly.SetEqConstraint(0, equalCons) ;
  _poly.SetEqConstant(0, 1.0) ;
  Vector currCons(initVariNum) ;
  for (int i = 0; i < colSize; ++ i) {
    currCons.head(initVariNum-1) = 
        - poly.get_coefficients().col(colIdx[i]).transpose() ;
    currCons(initVariNum-1) = 0.0 ;
    _poly.SetEqConstraint(i+1, currCons) ;
    _poly.SetEqConstant(i+1, 0.0) ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints tableau in float constructed" << std::endl ;
  _poly.Print() ;
  log_mtx.unlock() ;
#endif

  _constraint_num = _poly.get_eq_constraint_num() ;
  _variable_num = _poly.get_variable_num() ;

  _objective.f.setZero() ; 
  // Cx+b >= 0
  // construct lambda_1 to lambda_n
  for (int j = 0; j < _variable_num-1; ++ j) {
    _objective.f.col(j).head(_parameter_num) =
        - poly.get_coefficients().row(j) ;
    _objective.f(_parameter_num, j) = poly.GetConstant(j) ;
  }
  // construct lambda_0
  _objective.f(_parameter_num, initVariNum-1) = 1.0 ;
  
  // if all the coeffs of parameters are 0
  if ( _objective.f.block(0, 0, _parameter_num, initVariNum-1).isZero() ) {
    return ; 
  }
    
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "objective tableau in float constructed" << std::endl
      << _objective.f << std::endl ;
  log_mtx.unlock() ;
#endif

  RNumber curr ;
  // the first constraint is Cx = 1
  for (int j = 0; j < _variable_num; ++ j) {
    curr = GetRational( equalCons(j), cutScale ) ;
    _constraints.at(0, j) = curr ;
  }
  _constraints.at(0, colNum-1) = 1 ;
  for (int i = 0; i < _constraint_num-1; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      curr = GetRational(_poly.GetEqCoef(i+1, j), cutScale) ; 
      _constraints.at(i+1, j) = curr ; 
    }
    _constraints.at(i+1, colNum-1) = 0 ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints tableau in rational constructed" << std::endl ;
  print(_constraints) ;
  log_mtx.unlock() ;
#endif

  // add the first points to the point stack
  // The point [0,...,0,1] is just for test, it may be at the same with
  // central point, in which the algorithm will go wrong.
  Vector firstPoint(_parameter_num+1) ;
  RegionIdx nullRegion(-1, true) ;
  if (NB_INITIAL_POINTS == 0) {
    firstPoint.head(_parameter_num).setZero() ;
    firstPoint(_parameter_num) = 1.0 ;
    worklist.emplace_back(nullRegion, 0, std::move(firstPoint) ) ;
  } else {
    for(int k=0; k<NB_INITIAL_POINTS; k++) {
      do {
	firstPoint.head(_parameter_num).setRandom() ;
	firstPoint(_parameter_num) = 1.0 ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "first point: " << firstPoint << std::endl ;
  log_mtx.unlock() ;
#endif

      } while ( CheckFirstPoint(firstPoint.head(_parameter_num), _normalization_point) ) ;
      worklist.emplace_back(nullRegion, k, firstPoint ) ;
    }
  }

// choose the initial basis. Use the same form with GLPK
// i.e. make their coeffs as 1
// e.g. 3lambda1 + 2lambda2 + 4 lambda3 = 1 is equivalent to
// 3/4lambda1 + 2/4lambda2 + lambda3 = 1/4
// => 3/4lambda1 + 2/4lambda2 <= 1/4, and remove lambda3 in the objective function 
// glpk will do:
// lambda' = 3/4lambda1 + 1/2lambda2, lambda' <= 1/4, where lambda' = -lambda3 
// so the new matrix is 3/4 1/2 -1 1/4 
// TODO continue here
  bool feasible = ChooseBasis() ;
  if ( ! feasible) {
    _objective.f.setZero() ;
    _objective.r.set_zero() ;
  }

  //TODO continue here. Test it 
  bool notZero = MinimizeLp() ;
  if ( ! notZero) {
    _objective.f.setZero() ;
    _objective.r.set_zero() ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "plp constructor ends" << std::endl ;
  log_mtx.unlock() ;
#endif
}

// TODO
/*
void Plp::ProjSubEqs() {
  
}
*/

/*******************************************************************************
 * Construct plp for convex hull
*******************************************************************************/
Plp::Plp(Polyhedron& poly1,
         Polyhedron& poly2,
	 Worklist_t& worklist,
         optimal_container* optimals,
         const int NB_INITIAL_POINTS)
    : _constraints(poly1.get_variable_num()+2,
          poly1.get_constraint_num()+poly2.get_constraint_num()+3),
      _objective(poly1.get_variable_num()+1,
          poly1.get_constraint_num()+poly2.get_constraint_num()+3),
      _optimals(optimals), _optimals_count(0) {
  if ( poly1.get_variable_num() != poly2.get_variable_num() ) {
    std::cerr << "The two polyhedra are not in the same environment,"
        << " cannot compute the convex hull." << std::endl ;
    std::terminate() ;
  }

  _parameter_num = poly1.get_variable_num() ;
  // column num is _variable_num + 1, 
  // _parameter_num+1 as the last row is the constant
  int initVariNum = poly1.get_constraint_num()+poly2.get_constraint_num()+2 ;
  int colNum = initVariNum + 1 ;
  int initConsNum = poly1.get_variable_num()+2 ;
 
  // construct the constraints
  _constraints.set_zero() ;
  // to normalization, evaluate the objective a point 
  // which is inside the input polyhedron

  Point tmpPoint = GlpkInterface::GetCentralPoint(poly1) ;
  if ( tmpPoint.IsEmpty() ) {
    tmpPoint = GlpkInterface::GetSatPoint(poly1) ;
  }
  if( tmpPoint.IsEmpty() ) {
    std::cerr << "Minimization failed. Cannot find a point inside the polyhedra" << std::endl ;
    std::terminate() ;
  }

  int cutScale = GetScale( poly1.get_dis_lower_bound(), _parameter_num ) ;
  _normalization_point = GetSimplePoint( tmpPoint.get_coordinates(), cutScale) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Lower bound of diatance: " << poly1.get_dis_lower_bound() << std::endl ;
  log_mtx.unlock() ;
#endif

#ifdef VERIFY_PLP
  RMatrix rationalPoint(_parameter_num, 1) ;
  for (int i = 0; i < _parameter_num; ++ i) {
    rationalPoint.at(i, 0) = GetRational(_normalization_point(i), cutScale) ;
  }
  bool pointVerified = VerifyCentralPoint(rationalPoint, poly) ;
  assert(pointVerified) ;
#endif

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "tmp central point OK: " << tmpPoint.get_coordinates().transpose()
      << std::endl ;
  std::cout << "real central point OK: " << _normalization_point.transpose()
      << std::endl ;
  log_mtx.unlock() ;
#endif

  Vector equalCons(initVariNum) ;
  equalCons.setZero() ;
  
  int variNum1 = poly1.get_constraint_num() ;
  int variNum2 = poly2.get_constraint_num() ;
  equalCons.head(variNum1) =
      _normalization_point * -poly1.get_coefficients().transpose()
      + poly1.get_constants() ;

  equalCons(variNum1) = 1.0 ;
  _poly.SetSize(initConsNum, initVariNum) ;
  _poly.Init() ;
  _poly.SetConstraint( 0, equalCons) ;
  _poly.SetConstant(0, 1.0) ;
  //_poly.SetOperator(0, ConstraintOperator::equal) ;
  Vector currCons(initVariNum) ;

  // the matrix containts: lambda, lambda_0, lambda', lambda_0', constant
  for (int i = 0; i < poly1.get_variable_num(); ++ i) {
    currCons.head(variNum1) = 
        - poly1.get_coefficients().col(i).transpose() ;
    currCons(variNum1) = 0.0 ;
    currCons.segment(variNum1+1, variNum2) = 
        poly2.get_coefficients().col(i).transpose() ;
    currCons(initVariNum-1) = 0.0 ;
    _poly.SetConstraint(i+1, currCons) ;
    _poly.SetConstant(i+1, 0.0) ;
    //_poly.SetOperator(i+1, ConstraintOperator::equal) ;
  }
  currCons.head(variNum1) = poly1.get_constants() ;
  currCons(variNum1) = 1.0 ;
  currCons.segment(variNum1+1, variNum2) = -poly2.get_constants() ; 
  currCons(initVariNum-1) = -1.0 ;
  _poly.SetConstraint(initConsNum-1, currCons) ;
  _poly.SetConstant(initConsNum-1, 0.0) ;
  //_poly.SetOperator(initConsNum-1, ConstraintOperator::equal) ;
 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints tableau in float constructed" << std::endl ;
  _poly.Print() ;
  log_mtx.unlock() ;
#endif

  _constraint_num = _poly.get_constraint_num() ;
  _variable_num = _poly.get_variable_num() ;

  _objective.f.setZero() ; 
  // Cx+b >= 0
  // construct lambda_1 to lambda_n
  for (int j = 0; j < variNum1; ++ j) {
    _objective.f.col(j).head(_parameter_num) =
        - poly1.get_coefficients().row(j) ;
    _objective.f(_parameter_num, j) = poly1.GetConstant(j) ;
  }
  // construct lambda_0
  _objective.f(_parameter_num, variNum1) = 1.0 ;
    
  // if all the coeffs of parameters are 0
  if ( _objective.f.block(0, 0, _parameter_num, initVariNum-1).isZero() ) {
    return ; 
  }
    
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "objective tableau in float constructed" << std::endl
      << _objective.f << std::endl ;
  log_mtx.unlock() ;
#endif

  RNumber curr ;
  // the first constraint is Cx = 1
  for (int j = 0; j < _variable_num; ++ j) {
    curr = GetRational( equalCons(j), cutScale ) ;
    _constraints.at(0, j) = curr ;
  }
  _constraints.at(0, colNum-1) = 1 ;
  for (int i = 0; i < _constraint_num-1; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      curr = GetRational(_poly.GetCoef(i+1, j), cutScale) ; 
      _constraints.at(i+1, j) = curr ; 
    }
    _constraints.at(i+1, colNum-1) = 0 ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints tableau in rational constructed" << std::endl ;
  print(_constraints) ;
  log_mtx.unlock() ;
#endif

  // add the first points to the point stack
  // TODO remove this after testing. The point [0,...,0,1] may be at the same with
  // central point, in which the algorithm will go wrong.
  Vector firstPoint(_parameter_num+1) ;
  RegionIdx nullRegion(-1, true) ;
  if (NB_INITIAL_POINTS == 0) {
    firstPoint.head(_parameter_num).setZero();
    firstPoint(_parameter_num) = 1.0;
    worklist.emplace_back(nullRegion, 0, std::move(firstPoint) ) ;
  } else {
    for(int k=0; k<NB_INITIAL_POINTS; k++) {
      do {
	firstPoint.head(_parameter_num).setRandom() ;
	firstPoint(_parameter_num) = 1.0 ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "first point: " << firstPoint << std::endl ;
  log_mtx.unlock() ;
#endif

      } while ( CheckFirstPoint(firstPoint.head(_parameter_num), _normalization_point) ) ;
      //worklist.emplace( FromIndex(-1,k), firstPoint ) ;
      worklist.emplace_back(nullRegion, k, firstPoint ) ;
    }
  }

// choose the initial basis. Use the same form with GLPK
// i.e. make their coeffs as 1
// e.g. 3lambda1 + 2lambda2 + 4 lambda3 = 1 is equivalent to
// 3/4lambda1 + 2/4lambda2 + lambda3 = 1/4
// => 3/4lambda1 + 2/4lambda2 <= 1/4, and remove lambda3 in the objective function 
// glpk will do:
// lambda' = 3/4lambda1 + 1/2lambda2, lambda' <= 1/4, where lambda' = -lambda3 
// so the new matrix is 3/4 1/2 -1 1/4 
  bool feasible = ChooseBasis() ;
  if ( ! feasible) {
    _objective.f.setZero() ;
    _objective.r.set_zero() ;
  }

  bool notZero = MinimizeLp() ;
  if ( ! notZero) {
    _objective.f.setZero() ;
    _objective.r.set_zero() ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "plp constructor ends" << std::endl ;
  log_mtx.unlock() ;
#endif
}

/*******************************************************************************
 * This function will remove the variables whose coefficients are all 0
 * It may change the floating-point constraints _poly
*******************************************************************************/
/*
std::vector<int> Plp::PretreatConsF() {
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "PretreatConsF() starts." << std::endl ;
  log_mtx.unlock() ;
#endif
  std::vector<int> colIdx ;
  int consNum = _poly.get_eq_constraint_num() ;
  int variNum = _poly.get_variable_num() ;
  // keep lambda_0
  for (int j = 0; j < variNum-1; ++ j) {
    if ( ! _poly.get_coefficients().col(j).isZero() ) {
      colIdx.push_back(j) ;
    }
  }
  colIdx.push_back(variNum-1) ;
  int colSize = colIdx.size() ;
  // there are columns to be removed
  if (colSize != variNum) {
    Matrix oriMatrix( _poly.get_coefficients() ) ;
    Vector constants( _poly.get_constants() ) ;
    _poly.SetSize(consNum, colSize) ;
    _poly.Init() ;
    double curr ;
    for (int i = 0; i < consNum; ++ i) {
      for (int j = 0; j < colSize; ++ j) {
        curr = oriMatrix( i, colIdx[j] ) ;
        _poly.SetCoef(i, j, curr) ; 
      } 
    }
    _poly.set_constants(constants) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Constraints after removing columns:" << std::endl ;
  _poly.Print() ;
  log_mtx.unlock() ;
#endif
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "PretreatConsF() ends." << std::endl ;
  log_mtx.unlock() ;
#endif
  
  return colIdx ;
}
*/


/* Checks if the point is covered by one of the already seen regions.
 */

RegionIdx Plp::CheckCovered(const Point& point) {
#ifndef LOCKFREE_REGIONS
  std::lock_guard<std::mutex> guard(_optimals_mutex);
#endif

  if (_optimals_count > 0) {
    auto it = _optimals->begin();
    for (int i = 0; i < _optimals_count; ++i, ++ it) {
      int pos = it->region.f.PointOnBoundary(point) ;
      if ( pos == PointPos::inside ) {
        RegionIdx res(i) ;
        return res ;
      } 
      else if ( pos == PointPos::outside ) {
        continue ;
      }
      else {
        RegionIdx res(i, false, pos) ;
        return res ;
      } 
    }
  }
  RegionIdx res(-1, true) ;
  return res ;
}

/*******************************************************************************
 * Reconstructs the simplex tableau
 * @para constratins the tableau of constraints
 * @para obj the tableau of objective
 * @return the reconstructed tableau
*******************************************************************************/
RMatrix Plp::Reconstruct(const std::vector<int>& basicIdx) const {
  // avoid zero rows in cons
  if (_constraints.rows() != _constraint_num) {
    RMatrix newCons( Tool::GetBlock( _constraints, 0, 0, _constraint_num, _constraints.cols() ) ) ;
    return Reconstruct(basicIdx, newCons, _objective.r) ;
  }
  else {
    return Reconstruct(basicIdx, _constraints, _objective.r) ;
  }
}

RMatrix Plp::Reconstruct(const std::vector<int>& basicIdx, const RMatrix& cons, 
      const RMatrix& obj) const {
  int basicNum = basicIdx.size() ;
  RMatrix objMatrix(obj) ;
  int consNum = cons.rows() ;
  int variNum = cons.cols() - 1 ;

  std::vector<int> rowIdx ;
  bool isZeroRow ;
  for (int i = 0; i < consNum; ++ i) {
    isZeroRow = true ;
    for (int j = 0; j < basicNum; ++ j) {
      if( ! cons.at(i, basicIdx[j]).is_zero() ) {
        isZeroRow = false ;
        break ;
      }
    }
    if ( ! isZeroRow) {
      rowIdx.push_back(i) ;
    }
  }

  //RMatrix coeff(consNum, basicNum) ;
  RMatrix coeff(rowIdx.size(), basicNum) ;
  coeff.set_zero() ;
  RMatrix objCoeff(objMatrix.rows(), basicNum) ;
  objCoeff.set_zero() ;
  //for (int i = 0; i < consNum; ++ i) {
  for (unsigned i = 0; i < rowIdx.size(); ++ i) {
    for (int j = 0; j < basicNum; ++ j) {
      //coeff.at(i,j) = cons.at(i, basicIdx[j]) ;
      coeff.at(i,j) = cons.at(rowIdx[i], basicIdx[j]) ;
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

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "solve linear equation: " ;
  std::cout << " Lambda * " ; print(coeff) ;
  std::cout << " = " ; print(objCoeff) ;
  log_mtx.unlock() ;
#endif

  // dixon vs fraction_free
  //fraction-free seems to have better performance
  
  auto lambda( coeff.transpose().solve_fraction_free( objCoeff.transpose() ) ) ;
  // avoid zero rows in cons
 
  return RMatrix(objMatrix - lambda.transpose() * cons) ; 
}

/*******************************************************************************
 * remove duplicated columns and the columns only contain 0
*******************************************************************************/
std::vector<int> Plp::GetNonDupIdxByCol(const RMatrix& result,
    const std::vector<int>& idx) {
  bool dup ;
  std::vector<int> nonZeroIdx, polyIdx, satIdx ;
  int currIdx, cmpIdx ;
  for (unsigned i = 0; i < idx.size(); ++ i) {
    currIdx = idx[i] ;
    if ( ! Tool::IsZeroCol(result, currIdx) ) {
      nonZeroIdx.push_back(currIdx) ;
    }
  }
  if (nonZeroIdx.size() != 0) {
    polyIdx.push_back( nonZeroIdx[0] ) ;
    for (unsigned i = 1; i < nonZeroIdx.size(); ++ i) {
      dup = false ;
      currIdx = nonZeroIdx[i] ;
      for (int k = i-1; k >= 0; -- k) {
        cmpIdx = nonZeroIdx[k] ;
        if ( Tool::ColumnsEq(result, currIdx, cmpIdx) ) {
        
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Constraint " << currIdx << " is duplicated" << std::endl ;
  log_mtx.unlock() ;
#endif
  
          dup = true ;
          break ;
        }
      }
      if (! dup) {
        polyIdx.push_back(currIdx) ;
      }
    }
  } 
  return polyIdx ;
}

std::vector<int> Plp::GetNonDupIdxByRow(const RMatrix& result, 
      const std::vector<int>& idx) {
  RMatrix transposeM( result.transpose() ) ;
  return GetNonDupIdxByCol(transposeM, idx) ;
}

RMatrix Plp::GetSubMatrix(const RMatrix& ori, const std::vector<int>& idx) {
  unsigned int rows = ori.rows() ;
  RMatrix newMatrix( rows, idx.size() ) ;
  unsigned int curr ;
  for (unsigned int j = 0; j < idx.size(); ++ j) {
    curr = idx[j] ;
    for (unsigned int i = 0; i < rows; ++ i) {
      newMatrix.at(i, j) = ori.at(i, curr) ; 
    } 
  }
  return newMatrix ;
}

/*******************************************************************************
 * Extracts the result from the tableau and stores into _regions and _optimals
 * @para result the reconstructed tableau
 * @para glplInter the glpk interface which contains the information of simplex 
*******************************************************************************/
int Plp::StoreResult(const RMatrix& result, const std::vector<int>& basicIdx,
    const std::vector<int>& nonBasicIdx, bool testPara) {
  RMatrix currOptimal(_parameter_num+1, 1) ; 
  if (nonBasicIdx.size() == 0) {
    return RegionState::noCons ;
  }
  std::vector<int> nonMinRegionIdx = GetNonDupIdxByCol(result, nonBasicIdx) ;
  // zero region
  if (nonMinRegionIdx.size() == 0) {
    return RegionState::zero ;
  }
  RMatrix nonMinRegion = GetSubMatrix(result, nonMinRegionIdx) ; 
  double tmp ;
  Polyhedron currRegion(nonMinRegion.cols(), _parameter_num) ;
  for (int i = 0; i < nonMinRegion.cols(); ++ i) {
    int currIdx = nonMinRegionIdx[i] ;
    for (int j = 0; j < _parameter_num; ++ j) {
      tmp = -Tool::GetDouble( (RNumber)result.at(j, currIdx) ) ;
      currRegion.SetCoef( i, j, tmp) ;
    }
    tmp = Tool::GetDouble( (RNumber)result.at(_parameter_num, currIdx) ) ;
    currRegion.SetConstant( i, tmp) ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  print(result) ;
  for (int idx : nonBasicIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  std::cout << "curr region: " ; currRegion.Print() ;
  log_mtx.unlock() ;
#endif

  // check if the constraint is 1 >= 0
  bool isConstant = true ;
  for (int j = 0; j < _parameter_num; ++ j) {
    if ( ! result.at(j, result.cols()-1).is_zero() ) {
      isConstant = false ;
      break ;
    }
  }
  if ( ! isConstant) {
    for (int j = 0; j < _parameter_num+1; ++ j) {
      currOptimal.at(j, 0) = -result.at(j, result.cols()-1) ;
    }
  }
  
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  if ( ! isConstant) {
    std::cout << "curr optimal: " ; print(currOptimal) ;
  }
  else {
    std::cout << "curr optimal: 0 <= 1" << std::endl ;
  }
  log_mtx.unlock() ;
#endif

  Point interP = GlpkInterface::GetSatPoint(currRegion, 1) ;
  Polyhedron minimizedRegion ;
  std::vector<int> activeIdx ; 
  if ( interP.IsEmpty() ) {
    // if the region is flat
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "GetSatPoint() : Flat region, check in rational" << std::endl ;
  log_mtx.unlock() ;
#endif 

    std::vector<int> rowIdx ;
    for (int j = 0; j < result.rows(); ++ j) {
      if ( ! Tool::IsZeroRow(result, j) ) {
        rowIdx.push_back(j) ;
      }
    }
    RMatrix tmpRR(nonMinRegionIdx.size(), rowIdx.size()) ;
    for (unsigned i = 0; i < rowIdx.size(); ++ i) {
      for (unsigned j = 0; j < nonMinRegionIdx.size(); ++ j) {
        tmpRR.at(j,i) = result.at(rowIdx[i], nonMinRegionIdx[j]) ;
      }
    }

    bool isFlatR = IsFlatRegionR(tmpRR) ;
    if ( ! isFlatR ) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Float result is wrong. Region is not flat." << std::endl ;
  log_mtx.unlock() ;
#endif 

      activeIdx = MinimizePolyR(tmpRR) ;
      minimizedRegion = currRegion.GetSubPoly(activeIdx) ;
      //Don't set internal point and witness point
      //The checker at the end will deal with it
      minimizedRegion.set_internal_point( Vector() ) ;
      
    }
    else {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Float result is correct. Region is flat." << std::endl ;
  log_mtx.unlock() ;
#endif 

      return RegionState::flat ;
    }
  }
  else {

/*******************************tmp code********************************/

  // curr_ird_cons contains the index of the irredundant constraints in
  // the columns of result matrix

  if (testPara) {
    return RegionState::normal ;
  }
 
/*******************************tmp code ends***************************/



    currRegion.set_internal_point(interP) ;
    //minimizedRegion = currRegion.GetMinimizedPoly(true) ; 
    currRegion.Minimize(true, true) ;
    // checker starts
    if ( ! currRegion.unsureCons.empty() ) {
      // use rational solver
      std::vector<int> rowIdx ;
      for (int j = 0; j < result.rows(); ++ j) {
        if ( ! Tool::IsZeroRow(result, j) ) {
          rowIdx.push_back(j) ;
        }
      }
      RMatrix tmpRR(nonMinRegionIdx.size(), rowIdx.size()) ;
      for (unsigned i = 0; i < rowIdx.size(); ++ i) {
        for (unsigned j = 0; j < nonMinRegionIdx.size(); ++ j) {
          tmpRR.at(j,i) = result.at(rowIdx[i], nonMinRegionIdx[j]) ;
        }
      }

      for (unsigned i = 0; i < currRegion.unsureCons.size(); ++ i) {
        int currIdx = currRegion.unsureCons[i] ;
        bool redundant = ConsRedundantR(tmpRR, currIdx) ;
        if ( ! redundant ) {
          currRegion.Activate(currIdx) ;
          Vector interP = currRegion.get_internal_point().get_coordinates() ;
          Vector cons = currRegion.get_coefficients().row(currIdx) ;
          double distanceNum = interP * cons.transpose() + currRegion.GetConstant(currIdx) ;
          double distanceDen = cons.norm() ;
          double distance = distanceNum / distanceDen ;
          Vector witness = interP + cons * distance ;
          currRegion.AddWitnessPoint(currIdx, witness) ;
        }
      }
    }

    // the idx in activeIdx are sorted
    activeIdx = currRegion.GetActiveIdx() ;
    //minimizedRegion = currRegion.GetSubPoly(activeIdx) ;
    minimizedRegion = currRegion.GetMinimizedPoly(true, true) ; 

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "minimized region: " ; minimizedRegion.Print() ;
    log_mtx.unlock() ;
#endif 

  }

  int activeSize = activeIdx.size() ;
  RMatrix regionR(_parameter_num+1, activeSize) ;
  int actualIdx ;
  for (int i = 0; i < _parameter_num+1; ++ i) {
    for (int j = 0; j < activeSize; ++ j) {
      actualIdx = activeIdx[j] ;
      regionR.at(i, j) = nonMinRegion.at(i, actualIdx) ;  
    }
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "minimized region in rational: " << std::endl ;
  print(regionR.transpose()) ;
  log_mtx.unlock() ;
#endif

  Point internalP(minimizedRegion.get_internal_point());
  int optimalIndex;

  RMatrix* opt ;
  RMatrix zero(1,1) ;
  zero.set_zero() ;
  if ( ! isConstant ) {
    opt = &currOptimal ;
  }
  else {
    opt = &zero ;
  }

  std::vector<int> currIrdCons ;
  for (int idx : activeIdx) {
    currIrdCons.push_back( nonMinRegionIdx[idx] ) ;
  }

  Region reg( std::move(minimizedRegion), std::move(regionR), 
      std::move(internalP), std::move(basicIdx), std::move(currIrdCons) ) ;
  
  /* [DM] region enqueuing BEGIN */
  {
#if LOCKFREE_REGIONS
    // emplace_back of tbb vector returns the pointer to the new element,
    // so it is the same with distance from begin to end.
    optimalIndex = std::distance( _optimals->begin(),
        _optimals->emplace_back( std::move(reg), std::move(*opt) ) );
		    
#else
    std::lock_guard<std::mutex> guard(_optimals_mutex);
    optimalIndex = std::distance( _optimals->begin(), _optimals->end() );
    _optimals->emplace_back( std::move(reg), std::move(*opt) ) ;
    
#endif

#if LOCKFREE_REGIONS
    /* [DM]
       This loop updates the (completely enqueued) region count.
       It waits until the all previous regions have been completely enqueued.
    */
#if LOCKFREE_REGIONS_WAIT_FOR_CONDITION
    std::unique_lock<std::mutex> lck(_optimals_mutex);
#endif
    while (_optimals_count < optimalIndex) {
#if LOCKFREE_REGIONS_WAIT_FOR_CONDITION
      _optimals_cv.wait(lck);
#elif LOCKFREE_REGIONS_PAUSE
      _mm_pause();
#endif
      spin_count++;
    }
#endif
    assert(_optimals_count == optimalIndex);
    _optimals_count = optimalIndex+1;
#if LOCKFREE_REGIONS
#if LOCKFREE_REGIONS_WAIT_FOR_CONDITION
    _optimals_cv.notify_all();
#endif
#endif
  }
  /* [DM] region enqueuing END */

  // the constraints are: Ax + b >= 0
  // Note that the result of the optimal in the simplex tableau
  // is the inverse of the real optimal value
  // if the current point is not first point, check adjacency
  //RegionIdx indexInfo(optimalIndex) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "added region: " << optimalIndex << " into " ;
  if (isConstant) std::cout << "constant region" << std::endl ;
  else std::cout << "normal region" << std::endl ;
  log_mtx.unlock() ;
#endif

  return optimalIndex ;
}

void Plp::AddTasks(const Task& currTask, int optimalIndex, task_feeder& feeder) {
  RegionIdx indexInfo(optimalIndex) ;
  if (GetTaskPoint(currTask).size() != 0) {
    Point currPoint( GetTaskPoint(currTask).head(_parameter_num) ) ; 
    RegionIdx from = GetFromRegion(currTask) ;
    int frontier = GetFromFrontier(currTask) ;
    if ( ! IsNullRegion(from) ) {
      bool adjacent = AreAdjacent(from, frontier, indexInfo) ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    if (IsConstantRegion(from)) std::cout << "constant region " ;
    else std::cout << "normal region " ;
    std::cout << from.idx << " and " ;
    bool isConstant = IsConstantRegion(optimalIndex) ;
    if (isConstant) std::cout << "constant region " ;
    else std::cout << "normal region " ;
    std::cout << optimalIndex << " are adjacent? " << adjacent << std::endl ;
    log_mtx.unlock() ;
#endif

      if (! adjacent) {
        AddPoint(currPoint, indexInfo, from, frontier, feeder, true) ;
      }
    }
  }

  Polyhedron minimizedRegion = get_region_f(optimalIndex) ;

  // if the region is computed by rational solver
  if ( minimizedRegion.get_internal_point().IsEmpty() ) {
    return ;
  }

  std::vector<Point> witness(minimizedRegion.GetWitness());
  for (int i = 0; i < (int)witness.size(); ++ i) {
    Vector newPoint(_parameter_num+1) ;
    newPoint.head(_parameter_num) = witness[i].get_coordinates() ;
    newPoint(_parameter_num) = 1 ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "added witness point: " << newPoint << std::endl ;
  log_mtx.unlock() ;
#endif

    //feeder.add( Task( FromIndex(optimalIndex, i), std::move(newPoint) ) ) ;
    feeder.add( Task(indexInfo, i, std::move(newPoint) ) ) ;
  }
}

/*******************************************************************************
//TODO description
*******************************************************************************/
bool Plp::AreAdjacent(const RegionIdx& from, int frontier,
    const RegionIdx& currIdx) {
  const RMatrix& region1 = _optimals->at(from.idx).region.r ;
  const RMatrix& region2 = _optimals->at(currIdx.idx).region.r ;

  for (int j = 0; j < region2.cols(); ++ j) {
    bool res = Tool::ColumnsEq(region1, region2, frontier, j, true) ;
    if (res == true) { 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "AreAdjacent(): found constraint " << j << std::endl ;
  log_mtx.unlock() ;
#endif

      RMatrix pointR(_parameter_num, 1) ;
      pointR.set_zero() ;
      int varIdx = -1 ;
      RNumber val ; 
      for (int p = 0; p < _parameter_num; ++ p) {
        val = region1.at(p, frontier) ;
        if ( ! val.is_zero() ) {
          varIdx = p ;
          break ;
        }
      }
      if (varIdx == -1) {
        std::cerr << "AreAdjacent(): error, 0 constraint" << std::endl ;
        std::terminate() ;
      }
      pointR.at(varIdx, 0) = - region1.at(_parameter_num, frontier) / val ; 

      bool verify = CheckAdjacent(from, currIdx, pointR) ; 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "verify adjacency: " << verify << std::endl ;
  log_mtx.unlock() ;
#endif
      
      if (verify == true) {

/*******************************to be tested********************************/

#ifdef PLP_CHECKER_ON 
        set_region_adj(from.idx, frontier) ;
#endif

/*******************************to be tested ends***************************/

        return true ;
      }
    } 
  }
  return false ;
}

bool Plp::CheckAdjacent(const RegionIdx& from, const RegionIdx& currIdx,
    const RMatrix& point) const {
  RMatrix currOpt(1, _parameter_num), fromOpt(1, _parameter_num) ;
  RMatrix currTmp(1, 1), fromTmp(1, 1) ;
  RNumber currVal, fromVal ;
  bool fromIsConst = IsConstantRegion(from) ;
  bool currIsConst = IsConstantRegion(currIdx) ;
  if ( ! fromIsConst ) {
    const RMatrix& optRef = _optimals->at(from.idx).optimal ;
    for (int i = 0; i < _parameter_num; ++ i) {
      fromOpt.at(0, i) = optRef.at(i, 0) ;
    } 
    fromTmp = fromOpt * point ;
    fromVal = fromTmp.at(0, 0) + optRef.at(_parameter_num, 0) ;
  }
  else {
    fromVal.set_one() ;
  }
  if ( ! currIsConst ) {
    const RMatrix& optRef = _optimals->at(currIdx.idx).optimal ;
    for (int i = 0; i < _parameter_num; ++ i) {
      currOpt.at(0, i) = optRef.at(i, 0) ;  
    } 
    currTmp = currOpt * point ; 
    currVal = currTmp.at(0, 0) + optRef.at(_parameter_num, 0) ;
  }
  else {
    currVal.set_one() ; 
  }

  return fromVal == currVal ;
}

/*******************************************************************************
 * Add a new point if two regions are not adjacent 
*******************************************************************************/
void Plp::AddPoint(const Point& currPoint, const RegionIdx& curr,
    const RegionIdx& from, int fromConsIdx, 
    task_feeder& feeder, bool newRegion) const {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Add point between regions" << std::endl ;
  log_mtx.unlock() ;
#endif
 
  const Point& fromPoint = _optimals->at(from.idx).region.point ;

  if ( fromPoint.IsEmpty() ) return ;

  const Polyhedron& fromRegion = _optimals->at(from.idx).region.f ; 
  Ray rayOri(currPoint, fromPoint) ;
  double fromMiniDis = Raytracing::GetDistance(
        fromPoint,
        fromRegion.get_coefficients().row(fromConsIdx), 
        fromRegion.GetConstant(fromConsIdx),
        rayOri) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "curr point: " << currPoint << std::endl ;
  std::cout << "from point: " << fromPoint << std::endl ;
  log_mtx.unlock() ;
#endif
 
  double newDis ;
  double pointDis = currPoint.GetPointDistance(fromPoint) ; 
  int currConsIdx = -1 ;
  const Polyhedron& currRegion = _optimals->at(curr.idx).region.f ; 
  // for the points belong to known regions
  // if the point is not on a boundary
  int boundIdx = curr.on_boundary ;
  bool onBoundary = PointOnBoundary(curr) ;
  if ( ! onBoundary ) {
    Ray rayCurr(fromPoint, currPoint) ; 
    double currMiniDisInverse = 0 ;
    double currDisInverse ;
    for (int i = 0; i < (int)currRegion.get_constraint_num(); ++ i) {
      double eva = currRegion.GetConstant(i) -
          currRegion.get_coefficients().row(i) * 
          currPoint.get_coordinates().transpose() ;
      // for new region
      // test if on boundary
      if (newRegion) {
        double threshold = Tool::GetDotProductThreshold( currPoint.get_coordinates(),
        currRegion.get_coefficients().row(i) ) ;
        if ( Double::AreEqual(eva, 0, threshold) ) {
          onBoundary = true ;
          boundIdx = i ;
          break ;
        }
      }
      currDisInverse = Raytracing::GetDistanceInverse(
          currRegion.get_coefficients().row(i), 
          rayCurr, eva) ;
      if (currMiniDisInverse < currDisInverse) {
        currMiniDisInverse = currDisInverse ;
        currConsIdx = i ;
      }
    } 
    if ( ! onBoundary) {
      double currMiniDis ;
      if (currConsIdx == -1) {
        // new region and located on the boundary
        currMiniDis = 0 ;
      }
      else {
        currMiniDis = 1 / currMiniDisInverse ; 
      }

      // small region?
      // TODO change this
      double threshold = pointDis * std::numeric_limits<double>::epsilon() ;
      if (pointDis - fromMiniDis - currMiniDis < threshold) {

        // small region, or more likely the two regions are overlapping
        
#ifdef PRINT_WARNING
  log_mtx.lock() ;
        std::cout << "Warning: AddPoint(): cannot find position for new point."
            << " value: " << pointDis - fromMiniDis - currMiniDis << std::endl ;
  log_mtx.unlock() ;
#endif
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Warning: AddPoint(): cannot find position for new point."
      << " value: " << pointDis - fromMiniDis - currMiniDis << std::endl ;
  log_mtx.unlock() ;
#endif

        return ;
      }
      newDis = (pointDis + fromMiniDis - currMiniDis) / 2.0 ; 
    }
  }

  if (onBoundary) {
    newDis = (pointDis + fromMiniDis) / 2.0 ; 
    currConsIdx = boundIdx ;
  }
  Vector newPoint(_parameter_num+1) ;
  newPoint.head(_parameter_num) = fromPoint.get_coordinates() 
      + rayOri.get_direction() * newDis ;
  newPoint(_parameter_num) = 1 ;

  // check small region
  Vector checkP = newPoint.head(_parameter_num) ; 
  int pos = currRegion.PointOnBoundary(checkP) ;
  if (pos != PointPos::outside) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Warning: AddPoint(): new point"
      << newPoint << " is too close to the currRegion. Return." << std::endl ;
  log_mtx.unlock() ;
#endif
 
    return ;
  } 
  pos = fromRegion.PointOnBoundary(checkP) ;
  if (pos != PointPos::outside) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Warning: AddPoint(): new point"
      << newPoint << " is too close to the fromRegion. Return." << std::endl ;
  log_mtx.unlock() ;
#endif

    return ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "added new point: " << newPoint.transpose() << std::endl ;
  log_mtx.unlock() ;
#endif

  feeder.add( Task( from, fromConsIdx, std::move(newPoint) ) ) ;
}

/*******************************************************************************
 * Transfers double to rational number 
 * @para num the double number to be transfered
 * @return the rational number
*******************************************************************************/
RNumber Plp::GetRational(double num, int scale) const {
  if ( std::numeric_limits<long>::max() / pow(10, scale) < std::abs(num) ) {
    std::cerr << "GetRational(): float numberi " << num 
        << " is too large/small, cannot convert" << std::endl ;
    std::terminate() ;
  }
  double tmp = num * pow(10, scale) ;
  long int intNum = round(tmp) ; // some compilers don't know std::round
  long int intDen = pow(10, scale) ;
  long int gcd = Tool::GetGcd(intNum, intDen) ;
  ZNumber rationalNum(intNum / gcd) ;
  ZNumber rationalDen(intDen / gcd) ;
  RNumber rational(rationalNum, rationalDen) ;    
  return rational ;
}

/*******************************************************************************
 * Changes precision of the double number. The precision is the same with 
 * that of the function GetRational().  
 * @para num the double number to be cut
 * @para scale the scale of the number after cutting
 * @return the double number with lower precision
*******************************************************************************/
double Plp::CutDouble(double num, int scale) const {
  if ( std::numeric_limits<double>::max() / pow(10, scale) < std::abs(num) ) {
    std::cerr << "CutDouble(): float number " << num
        << " is too large/small, cannot convert" << std::endl ;
    std::terminate() ;
  }
  double tmp = num * pow(10, scale) ;
  double newDouble = round(tmp) ; // some compilers don't know std::round
  return newDouble / pow(10, scale) ;
}

/*******************************************************************************
 * Minimize polyhedron with plp  
 * @return the minimized polyhedron
 * TODO maybe remove it later
*******************************************************************************/

/*******************************************************************************
 * check if the first point is the same with one of the normalization point 
 * (i.e. if the point will make the LP objective will be the same with one of
 * its constraints)
 * We need to ignore the variables which are projected
 * (They should NOT be the same)
 * @para point the point to be checked
 * @return true if the two points are the same
*******************************************************************************/
bool Plp::CheckFirstPoint(const Vector& central_point,
    const Vector& first_point) const {
  for (int i = 0; i < central_point.size(); ++ i) {
    if (std::find(_proj_idx.begin(), _proj_idx.end(), i) != _proj_idx.end()) {
      continue ;
    }
    else if (central_point[i] != first_point[i]) {
      return false ;
    }
  }
  return true ; 
}

int Plp::GetScale(double bound, int dimension) {
  int cutScale = 1 ;
  // if disLowerBound >= 1 then cutScale = 1.
  while (bound < 1.0) {
    bound *= pow(10, cutScale) ;
    cutScale += 1 ;
  }
  // when the precision is larger than 100
  if (dimension > 100) {
    cutScale += 1 ;
  }
  return cutScale ;
}

Vector Plp::GetSimplePoint(const Vector& point, int scale) {
  Vector realPoint( point.size() ) ;
  for (int i = 0; i < point.size(); ++ i) {
    realPoint(i) = CutDouble(point(i), scale) ;
  }
  return realPoint ;
}

// choose the last n columns as the initial basic variables
// TODO this function is a mess, change it
/*******************************************************************************
 * @results initialize _constraint_num, _variable_num and the matrix
 * _constraints and matrix _objective
 * @return true if the constraints are feasible
*******************************************************************************/
bool Plp::ChooseBasis() {
  const int initVariNum = _constraints.cols() - 1 ;
  const int constraintsCols = _constraints.cols() ;

  ReducedMatrix reducedM = Tool::GaussianElimination(_constraints) ;
  std::vector<int> rowOrder = reducedM.rowIdx ;
  int rowSize = rowOrder.size() ;
  std::vector<int> basisOrder = reducedM.colIdx ;
  RMatrix origen( std::move(reducedM.matrix) ) ;

  //put the basic variables at last for constraints in rational
  //RMatrix origen(_constraints) ;
  _constraints.set_zero() ;
  for (int i = 0; i < rowSize; ++ i) {
    for (int j = 0; j < initVariNum; ++ j) {
      _constraints.at(i, j) = origen.at( rowOrder[i], basisOrder[j] ) ;
    }
    _constraints.at(i, initVariNum) = origen.at(rowOrder[i], initVariNum) ;
  } 

  _constraint_num = rowSize ;
  // _variable_num does not change
  int nonBasicNum = _variable_num - _constraint_num ;
  // now the top-left block of _constraints has no zero rows
  // and we obtained the matrix in echelon form 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "constraints after moving columns: " << std::endl ; 
  print(_constraints) ;
  log_mtx.unlock() ;
#endif

  // save the order of columns of objective
  std::vector<int> objColOrder = basisOrder ;
 
  bool initConsFeasible = true ;
  RNumber zero ;
  for (int i = 0; i < _constraint_num; ++ i) {
    if (_constraints.at(i, constraintsCols-1) < zero) {
      initConsFeasible = false ;
      break ;
    }
  }

  if ( ! initConsFeasible ) {
    Simplex solver(_constraints, _constraint_num, _variable_num) ;  
    std::vector<int> feasibleIdx = solver.Solve(true) ;
    // if the problem is infeasible
    if ( feasibleIdx.empty() ) {
      return false ;
    }

    std::vector<int> feasibleAll ;
    for (int j = 0; j < _variable_num; ++ j) {
      if ( std::find(feasibleIdx.begin(), feasibleIdx.end(), j) == feasibleIdx.end() ) {
        feasibleAll.push_back(j) ;
      }
    }
    for (unsigned j = 0; j < feasibleIdx.size(); ++ j) {
      feasibleAll.push_back( feasibleIdx[j] ) ;
    }
    RMatrix feasibleCons( solver.get_cons_feasible() ) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "feasible order: " ;
  for (unsigned j = 0; j < feasibleAll.size(); ++ j) {
    std::cout << " " << feasibleAll[j] << " " ;
  }
  std::cout << std::endl ;
  std::cout << "feasible constraints: "; print(feasibleCons) ;
  log_mtx.unlock() ;
#endif
 
    //put the basic variables at last for constraints in rational
    _constraints.set_zero() ;
    for (int i = 0; i < _constraint_num; ++ i) {
      for (int j = 0; j < _variable_num; ++ j) {
        _constraints.at(i, j) = feasibleCons.at( i, feasibleAll[j] ) ;
      }
      _constraints.at(i, initVariNum) = feasibleCons.at(i, initVariNum) ;
    } 

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "feasible constraints after moving columns: " << std::endl ; 
    print(_constraints) ;
    log_mtx.unlock() ;
#endif

    // change the columns order of objective
    std::vector<int> tmpOrder = objColOrder ;
    assert( feasibleAll.size() == objColOrder.size() ) ;
    for (unsigned i = 0; i < tmpOrder.size(); ++ i) {
      objColOrder[i] = tmpOrder[ feasibleAll[i] ] ;
    }
  }
#ifdef DEBUGINFO_PLP
  else {
    log_mtx.lock() ;
    std::cout << "Initial basis is feasible." << std::endl ; 
    log_mtx.unlock() ;
  }
#endif

  // change the order of objective in float
  Matrix feasibleObj = _objective.f ;
  for (int j = 0; j < _variable_num; ++ j) {
    _objective.f.col(j) = feasibleObj.col( objColOrder[j] ) ;
  } 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Objective in float: "
      << std::endl << _objective.f << std::endl ;
  log_mtx.unlock() ;
#endif
  
  _objective.r.set_zero() ;
  long int val ;
  for (int i = 0; i < _parameter_num+1; ++ i) {
    for(int j = 0; j < constraintsCols; ++ j) { 
      val = Tool::GetIntOfFloat( _objective.f(i, j) ) ;
      _objective.r.at(i, j) = val ;
    }
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Changed variables order of objective in rational: "
      << std::endl ;
  print(_objective.r) ;
  log_mtx.unlock() ;
#endif

  RMatrix coeff(_constraint_num, _constraint_num) ;
  coeff.set_zero() ;
  RMatrix objCoeff(_parameter_num+1, _constraint_num) ;
  objCoeff.set_zero() ;
  
  for (int i = 0; i < _constraint_num; ++ i) {
    coeff.at(i,i) = 1 ;
  }

  for (int i = 0; i < _parameter_num+1; ++ i) {
    for (int j = 0; j < _constraint_num; ++ j) {
      objCoeff.at(i,j) = _objective.r.at(i, nonBasicNum+j) ;
    }
  }

  auto lambda( coeff.transpose().solve_fraction_free( objCoeff.transpose() ) ) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "solve linear equation: " ;
  std::cout << " Lambda * " ; print(coeff) ;
  std::cout << " = " ; print(objCoeff) ;
  log_mtx.unlock() ;
#endif

  if (lambda.rows() != _constraints.rows()) {
    RMatrix newCons( Tool::GetBlock(_constraints, 0, 0,
        lambda.rows(), constraintsCols) ) ;
    _objective.r = _objective.r - lambda.transpose() * newCons ; 
  }
  else {
    _objective.r = _objective.r - lambda.transpose() * _constraints ; 
  }

  double currF ;
  _objective.f.setZero() ;
  for (int i = 0; i < _parameter_num+1; ++ i) {
    for (int j = 0; j < constraintsCols; ++ j) {
      currF = Tool::GetDouble( (RNumber)_objective.r.at(i, j) ) ;
      _objective.f(i, j) = currF ;
    }
  } 

  // remove the constraints which only contain the constants, i.e. lambda_i = c
  std::vector<int> ineqIdx ;
  for (int i = 0; i < _constraint_num; ++ i) {
    if ( ! Tool::IsZeroBlock(_constraints, i, 1, 0, nonBasicNum) ) {
      ineqIdx.push_back(i) ;
    }
  }
  if ( (int)ineqIdx.size() != _constraint_num ) {
    int removed = _constraint_num - ineqIdx.size() ;
    _constraint_num = ineqIdx.size() ;
    // the number of basic variable = the number of constraints
    _variable_num -= removed ;
    // nonBasicNum does not change
    RMatrix oriCons(_constraints) ;
    _constraints.set_zero() ;
    for (int i = 0; i < _constraint_num; ++ i) {
      for (int j = 0; j < nonBasicNum; ++ j) {
        _constraints.at(i, j) = oriCons.at(ineqIdx[i], j) ;
      } 
      _constraints.at(i, i+nonBasicNum) = 1 ;
      _constraints.at(i, initVariNum) = oriCons.at(ineqIdx[i], initVariNum) ;
    }
  }

  if (_poly.get_constraint_num() != _constraint_num) {
    _poly.SetSize(_constraint_num, _variable_num) ;
    _poly.Init() ;
  }
  for (int i = 0; i < _constraint_num; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      currF = Tool::GetDouble( (RNumber)_constraints.at(i, j) ) ;
      _poly.SetCoef(i, j, currF) ;
    }
    currF = Tool::GetDouble( (RNumber)_constraints.at(i, constraintsCols-1) ) ;
    _poly.SetConstant(i, currF) ;
    //_poly.SetOperator(i, ConstraintOperator::equal) ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Feasible constraints in float: " ;
  _poly.Print() ;
  log_mtx.unlock() ;
#endif

  _glpk_poly.SetSize(_constraint_num, nonBasicNum) ;
  _glpk_poly.set_coefficients( _poly.get_coefficients().leftCols(nonBasicNum) ) ;
  _glpk_poly.set_constants( _poly.get_constants() ) ;

  for (int j = 0; j < nonBasicNum; ++ j) {
    _non_negative_idx.push_back(j) ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "glpk poly: " << std::endl ; _glpk_poly.Print() ;
  log_mtx.unlock() ;
#endif

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "initial obj in rational: " << std::endl ; print(_objective.r) ;
  std::cout << "initial obj in float: " << std::endl << _objective.f << std::endl ;
  std::cout << "final constraints in rational: " << std::endl ;
  print(_constraints) ;
  std::cout << "cons_num: " << _constraint_num << " vari_num: "
      << _variable_num << std::endl ;
  log_mtx.unlock() ;
#endif

  return true ;
}

/*******************************************************************************
 * This function is only used in the function MinimizedLp().
 * @return the index of the variable x if the row in the matrix is in the form 
 * of -ax <= 0, where a is a positive constant; otherwise return -1.
*******************************************************************************/
int Plp::GetBdIdx(const RMatrix& matrix, int rowIdx) {
  unsigned colNum = matrix.cols() ;
  RNumber zero ;
  zero.set_zero() ;
  bool isBd = false ;
  int bdIdx ;
  for (unsigned j = 0; j < colNum; ++ j) {
    if (matrix.at(rowIdx,j) > zero) {
      bdIdx = -1 ;
      break ;
    }
    if (matrix.at(rowIdx,j) < zero) {
      if ( ! isBd) {
        bdIdx = j ;
        isBd = true ;
      }
      else {
        isBd = false ;
        bdIdx = -1 ;
        break ;
      }
    }
  }
  return bdIdx ;
}


/*******************************************************************************
 * Minimize the constraints of the LP problem (the polyhedron of lambda)
 * We need to add the constraints lambda_i >= 0
 * Before raytracing, we need to remove the duplicated constraints in ratioanl
 * constraints, and then convert the non-duplicated constraints into floating
 * point numbers
 * There are two steps: first we remove the implicit equations, and then
 * minimize the polyhedron which does not contain the equalities.
 * @results the minimized polyhedron saved in _constraints in rational
 * and _glpk_poly in floating-point
 * @return true if the constraints construct a nonempty-interior poly
*******************************************************************************/
bool Plp::MinimizeLp() {
  // if the initial problem has no non-basic variables
  if (_variable_num == _constraint_num) return true ;
  int rowNum = _variable_num ;
  int colNum = _variable_num - _constraint_num + 1 ;
  const int constraintsCols = _constraints.cols() ;

  RMatrix tmpConstraints(rowNum, colNum) ;
  // copy coefficients
  Tool::CopyMatrix(_constraints, tmpConstraints, 0, 0, 0, 0,
      _constraint_num, colNum-1) ; 
  // copy constants
  Tool::CopyMatrix(_constraints, tmpConstraints, 0, constraintsCols-1,
      0, colNum-1, _constraint_num, 1) ;
  // add the constraints lambda_i >= 0
  for (int i = 0; i < colNum-1; ++ i) {
    tmpConstraints.at(_constraint_num+i, i) = -1 ;
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "The LP matrix to be minimized:" << std::endl ;
  print(tmpConstraints) ;
  log_mtx.unlock() ;
#endif

  std::vector<int> idx ;
  for (int i = 0; i < rowNum; ++ i) {
    idx.push_back(i) ; 
  }
  std::vector<int> nonDupIdx = GetNonDupIdxByRow(tmpConstraints, idx) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "nonDupIdx: " ;
  for (int i : nonDupIdx) {
    std::cout << i << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

  std::vector<int> currActiveIdx = nonDupIdx ;
  bool toMinimize = true ;
  Point point ;
  Polyhedron polyToMini ;
  
  std::vector<int> removedVari ;
  // loop until there is no implicit constraints
  while (true) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "The loop of searching equalities starts" << std::endl ;
  log_mtx.unlock() ;
#endif

    Polyhedron tmpPoly(tmpConstraints, currActiveIdx) ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "poly: " ; tmpPoly.Print() ;
    log_mtx.unlock() ;
#endif

    point = GlpkInterface::GetCentralPoint(tmpPoly, true) ;
    if ( ! tmpPoly.IsFlat() ) {

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "The poly is not flat" << std::endl ;
    log_mtx.unlock() ;
#endif

      polyToMini = std::move(tmpPoly) ; 
      break ;
    }

    // if the poly is flat, start to find the implicit constraints 
    
    // cannot find a point
    // the poly may be infeasible, but we cannot trust floating-point computation
    // just give up minimization and continue, i.e. return true.
    if ( point.IsEmpty() ) {

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "MinimizeLp() cannot find a interior point" << std::endl ;
  log_mtx.unlock() ;
#endif

      return true ;
    }

    // if there are implicit equations
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "MinimizeLp(): the poly is flat, "
      << "start to find implicit equations with the point: " << point << std::endl ;
  log_mtx.unlock() ;
#endif

    // get the constraints which construct the implicit equation
    std::vector<int> implicitIdx ;
    double tmp ;
    for (int i  = 0; i < tmpPoly.get_constraint_num(); ++ i) {
      tmp = tmpPoly.get_coefficients().row(i) * point.get_coordinates().transpose()
          - tmpPoly.GetConstant(i) ;
      if ( Double::AreEqual(tmp, 0, 1e-6) ) {
        implicitIdx.push_back( currActiveIdx[i] ) ;
      }
    }

    // it should be empty
    if ( implicitIdx.empty() ) {
      std::cerr << "MinimizeLp(): cannot find implicit constraints" << std::endl ;
      std::terminate() ;
    }

    int impIdxSize = implicitIdx.size() ;
    RMatrix calFarkas(colNum, impIdxSize) ;
    std::vector<int> basicIdx, oriIdx ;
    std::vector<RNumber> feasibleV ;
    std::vector<int> eqIdx ;
    unsigned eqConsNum, consColSize ;
    int currIdx ;
    bool found = false ;
    // this loop is used to avoid the wrong constraints found by floating-point
    // in most cases we can get the combination when k=0
    for (int k = 0; k < impIdxSize; ++ k) {
      feasibleV.clear() ;
      eqIdx.clear() ;
      basicIdx.clear() ;
      oriIdx.clear() ;
      // if the kth constraint is the (negative) combination of others
      for (int i = 0, m = 0; i < impIdxSize; ++ i) {
        if (i == k) continue ;
        for (int j = 0; j < colNum-1; ++ j) {
          calFarkas.at(j,m) = tmpConstraints.at(implicitIdx[i],j) ;
        }
        calFarkas.at(colNum-1,m) = - tmpConstraints.at(implicitIdx[i], colNum-1) ;
        ++ m ;
      } 
      for (int j = 0; j < colNum; ++ j) {
        calFarkas.at(j,impIdxSize-1) = - tmpConstraints.at(implicitIdx[k], j) ;
      }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "The constraints of Farkas: " << std::endl ;
    print(calFarkas) ;
    log_mtx.unlock() ;
#endif

      // call Gaussian for echelon form
      // TODO move the Gaussian into Simplex solver
      ReducedMatrix echelonM = Tool::GaussianElimination(calFarkas) ;
      // if the problem is infeasible, e.g. 0 = -2
      if ( ! echelonM.succeed) {
        break ;
      }
      std::vector<int> rowOrder = echelonM.rowIdx ;
      int rowSize = rowOrder.size() ;
      int colSize = echelonM.colIdx.size() ;
      int initColSize = echelonM.matrix.cols() ;
      std::vector<int> basisOrder = echelonM.colIdx ;
      RMatrix reducedM(rowSize, colSize+1) ;
      for (int i = 0; i < rowSize; ++ i) {
        for (int j = 0; j < colSize; ++ j) {
          reducedM.at(i, j) = echelonM.matrix.at( rowOrder[i], basisOrder[j] ) ;
        }
        reducedM.at(i, colSize) = echelonM.matrix.at(rowOrder[i], initColSize-1) ;
      }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "The poly in echelon form: " << std::endl ;
    print(reducedM) ;
    log_mtx.unlock() ;
#endif

      Simplex solver(reducedM, rowSize, colSize) ;  
      std::vector<int> feasibleIdx = solver.Solve(true) ;
      // if the problem is infeasible, i.e. the current constraint doesn't construct the combination
      // this is caused by the inaccuracy of floating point
      if ( feasibleIdx.empty() ) {
        continue ;
      }
      currIdx = k ;
      RMatrix feaV( solver.get_feasible_vertex() ) ;

      // the constraints construct the combination 
      for (unsigned i = 0; i < feasibleIdx.size(); ++ i) {
        if ( ! feaV.at(0,i).is_zero() ) {
          basicIdx.push_back( basisOrder[ feasibleIdx[i] ] ) ;
          feasibleV.push_back( (RNumber)feaV.at(0, i) ) ;
        }
      }

      // all the constraints in the farkas except k
      for (int i = 0; i < impIdxSize; ++ i) {
        if (i != k) {
          oriIdx.push_back( implicitIdx[i] ) ;
        }
      }

      for (int idx : basicIdx) {
        eqIdx.push_back( oriIdx[idx] ) ;
      }
      eqIdx.push_back(implicitIdx[currIdx]) ;
      RNumber one ;
      one.set_one() ;
      feasibleV.push_back(one) ;
      eqConsNum = feasibleV.size() ;
      consColSize = tmpConstraints.cols() ;
      RMatrix eqMatrix(eqConsNum, consColSize) ;
      for (unsigned i = 0; i < eqConsNum; ++ i) {
        for (unsigned j = 0; j < consColSize; ++ j) {
          eqMatrix.at(i,j) = tmpConstraints.at(eqIdx[i],j) ;
        }
      }
      ReducedMatrix reducedEqM = Tool::GaussianElimination(eqMatrix) ;

      // until here we found the combination
      if (reducedEqM.succeed) {
        found = true ;
        break ;
      }

    }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "finding combination ends." << std::endl ;
    log_mtx.unlock() ;
#endif

    // if the poly is still flat after removing the implicit constraints
    if ( ! found) {
      toMinimize = false ;
      break ;
    }

    // we found a set of inequalities which construct thei implicit equalities
    // now we start to extract the equalities
    
    // verification
#ifdef VERIFY_PLP
    RMatrix verify(1, consColSize) ;
    verify.set_zero() ;
    for (unsigned i = 0; i < eqConsNum; ++ i) {
      for (unsigned j = 0; j < eqConsNum; ++ j) {
        verify.at(0,j) += feasibleV[i] * tmpConstraints.at(eqIdx[i],j) ;
      }
    }
    assert( verify.is_zero() ) ;
#endif
    // verification ends
    
    RMatrix eqMatrix(eqConsNum, consColSize) ;
    for (unsigned i = 0; i < eqConsNum; ++ i) {
      for (unsigned j = 0; j < consColSize; ++ j) {
        eqMatrix.at(i,j) = tmpConstraints.at(eqIdx[i],j) ;
      }
    }
    ReducedMatrix reducedEqM = Tool::GaussianElimination(eqMatrix) ;
    std::vector<int> eqRowOrder = reducedEqM.rowIdx ;
    std::vector<int> eqBasicIdx = reducedEqM.basicIdx ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
    std::cout << "Matrix of equalities" << std::endl ;
    print(eqMatrix) ;
    std::cout << "Reduced matrix of equalities" << std::endl ;
    print(reducedEqM.matrix) ;
  log_mtx.unlock() ;
#endif

    unsigned ineqConsSize = tmpConstraints.rows() - eqConsNum ; 
    RMatrix ineqMatrix(ineqConsSize, consColSize) ;
    std::vector<int> ineqIdx ;
    for (int i = 0; i < tmpConstraints.rows(); ++ i) {
      if ( std::find(eqIdx.begin(), eqIdx.end(), i) == eqIdx.end() ) {
        ineqIdx.push_back(i) ;
      }
    }
    assert(ineqIdx.size() == ineqConsSize) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "The index of inequalities in the matrix of poly to be minimized: " ;
  for (int idx : ineqIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

    for (unsigned i = 0; i < ineqConsSize; ++ i) {
      for (unsigned j = 0; j < consColSize; ++ j) {
        ineqMatrix.at(i,j) = tmpConstraints.at(ineqIdx[i],j) ;
      }
    }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
    std::cout << "Matrix of inequalities" << std::endl ;
    print(ineqMatrix) ;
  log_mtx.unlock() ;
#endif
  
    for (unsigned i = 0; i < eqRowOrder.size(); ++ i) {
      RMatrix ratio( Tool::GetBlock(ineqMatrix, 0, eqBasicIdx[i], ineqConsSize, 1) ) ;
      RMatrix currEq(Tool::GetBlock(reducedEqM.matrix, eqRowOrder[i], 0, 1, consColSize) ) ;
      ineqMatrix -= ratio * currEq ; 
      removedVari.push_back( eqBasicIdx[i] ) ;
    }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
    std::cout << "Reduced matrix of inequalities" << std::endl ;
    print(ineqMatrix) ;
  log_mtx.unlock() ;
#endif
  
    bool isZero = Tool::IsZeroBlock(ineqMatrix, 0, ineqConsSize, 0, consColSize-1) ;
    if (isZero) {
      return false ;
    }

    // evaluate the variable in the objective
    int objRowSize = _objective.r.rows() ;
    int objColSize = _objective.r.cols() ;
    RMatrix tmpEqM(1, objColSize) ;
    tmpEqM.set_zero() ;
    for (unsigned i = 0; i < eqRowOrder.size(); ++ i) {
      RMatrix ratio( Tool::GetBlock(_objective.r, 0, eqBasicIdx[i], objRowSize, 1) ) ;
      for (unsigned j = 0; j < consColSize-1; ++ j) {
        tmpEqM.at(0,j) = reducedEqM.matrix.at(eqRowOrder[i],j) ;
      }
      tmpEqM.at(0, objColSize-1) = reducedEqM.matrix.at(eqRowOrder[i], consColSize-1) ;
      
      _objective.r -= ratio * tmpEqM ; 
    }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
    std::cout << "Reduced objective" << std::endl ;
    print(_objective.r) ;
  log_mtx.unlock() ;
#endif

    std::vector<int> tmpIdx ;
    for (unsigned i = 0; i < ineqMatrix.rows(); ++ i) {
      tmpIdx.push_back(i) ;
    } 
    std::vector<int> nonDupIdx = GetNonDupIdxByRow(ineqMatrix, tmpIdx) ;
    std::vector<int> ineqIdxRef ;
    for (unsigned i = 0; i < nonDupIdx.size(); ++ i) {
      ineqIdxRef.push_back( ineqIdx[ nonDupIdx[i] ] ) ;
    }
    tmpConstraints.set_zero() ;
    for (unsigned i = 0; i < nonDupIdx.size(); ++ i) {
      for (unsigned j = 0; j < consColSize; ++ j) {
        tmpConstraints.at(ineqIdxRef[i],j) = ineqMatrix.at(nonDupIdx[i],j) ;
      }
    }
    currActiveIdx = std::move(ineqIdxRef) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "The reduced matrix of the poly to be minimized: "  << std::endl ;
  print(tmpConstraints) ;
  std::cout << "The index of  remained constraints: " ;
  for (int idx : currActiveIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  std::cout << "The loop of searching equalities ends" << std::endl ;
  log_mtx.unlock() ;
#endif

  }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "currActiveIdx: " ;
    for (int i : currActiveIdx) {
      std::cout << i << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif

  if (toMinimize) {
    polyToMini.set_internal_point(point) ;
    polyToMini.Minimize() ;
    if ( ! polyToMini.unsureCons.empty() ) {
      // there is the constraint 0 <= 1
      unsigned activeSize = currActiveIdx.size() ;
      unsigned consCols = tmpConstraints.cols() ;
      RMatrix tmpRR(activeSize+1, consCols) ;
      for (unsigned i = 0; i < activeSize; ++ i) {
        for (unsigned j = 0; j < consCols; ++ j) {
          tmpRR.at(i,j) = tmpConstraints.at(currActiveIdx[i], j) ;
        }
      }
      tmpRR.at(activeSize, consCols-1).set_one() ;
      for (unsigned i = 0; i < polyToMini.unsureCons.size(); ++ i) {
      int currIdx = polyToMini.unsureCons[i] ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "Check unsure: " << currIdx << std::endl ;
    log_mtx.unlock() ;
#endif

        bool redundant = ConsRedundantR(tmpRR, currIdx) ;
        if ( ! redundant ) {
          polyToMini.Activate(currIdx) ;

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "irredundant" << std::endl ;
    log_mtx.unlock() ;
#endif

        }
      }

    }

    std::vector<int> minimizedIdx = polyToMini.GetActiveIdx() ;
    // the value in activeIdx is the index of index in nonDupIdx
    std::vector<int> tmpIdx = std::move(currActiveIdx) ;
    currActiveIdx.clear() ;
    for (unsigned i = 0; i < minimizedIdx.size(); ++ i) {
      currActiveIdx.push_back( tmpIdx[minimizedIdx[i]] ) ;
    }
    

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "minimizedIdx: " ;
    for (int i : minimizedIdx) {
      std::cout << i << " " ;
    }
    std::cout << std::endl ;
    std::cout << "currActiveIdx: " ;
    for (int i : currActiveIdx) {
      std::cout << i << " " ;
    }
    std::cout << std::endl ;
    log_mtx.unlock() ;
#endif
  }

  int activeSize = currActiveIdx.size() ;
  // if no redundant constraints
  if (activeSize == rowNum) {
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Do not need to be minimized" << std::endl ;
  log_mtx.unlock() ;
#endif

    return true ;
  }
  std::vector<int> consIdx, bdIdx ;
  int currIdx ;
  int variIdx ;
  for (int i = 0; i < activeSize; ++ i) {
    currIdx = currActiveIdx[i] ;
    variIdx = GetBdIdx(tmpConstraints, currIdx) ;
    if (variIdx == -1) {
      consIdx.push_back(currIdx) ;
    }
    else {
      bdIdx.push_back(variIdx) ;
    }
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "ConsIdx: " ;
  for (int idx : consIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  std::cout << "bdIdx: " ;
  for (int idx : bdIdx) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

  int newConsNum = consIdx.size() ;
  int newBdNum = bdIdx.size() ;
  std::vector<int> remainedVari ;
  for (unsigned i = 0; i < tmpConstraints.cols()-1; ++ i) {
    if ( std::find(removedVari.begin(), removedVari.end(), i) == removedVari.end() ) {
      remainedVari.push_back(i) ;
    }
  }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Removed variables: " ;
  for (int idx : removedVari) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  std::cout << "remainedVari: " ;
  for (int idx : remainedVari) {
    std::cout << idx << " " ;
  }
  std::cout << std::endl ;
  log_mtx.unlock() ;
#endif

  if ( newConsNum < _constraint_num || ! removedVari.empty() ) {
    // the variables removed by substituting one variable in the equalities
    _variable_num -= removedVari.size() ;
    // When we remove a constraint, one basic variable is also removed
    int removed = _constraint_num - newConsNum ;
    _variable_num -= removed ;
    _constraint_num = newConsNum ;
    std::sort( consIdx.begin(), consIdx.end() ) ;
    //Put the non-zero sub-matrix at the top-left
    _constraints.set_zero() ;
    int nonBasicNum = _variable_num - _constraint_num ;
    assert( nonBasicNum == (int)remainedVari.size() ) ;
    // non-basic variables
    for (int i = 0; i < newConsNum; ++ i) {
      for (int j = 0; j < nonBasicNum; ++ j) {
        _constraints.at(i,j) = tmpConstraints.at(consIdx[i],remainedVari[j]) ; 
      }
    }
    // basic variables
    for (int i = 0 ; i < newConsNum; ++ i) {
      _constraints.at(i, i+nonBasicNum) = 1 ;
    }
    // constants
    for (int i = 0; i < newConsNum; ++ i) {
      _constraints.at(i, constraintsCols-1)
          = tmpConstraints.at(consIdx[i], tmpConstraints.cols()-1) ;
    } 

    // constraints in floating point
    _glpk_poly.SetSize(_constraint_num, nonBasicNum) ;
    _glpk_poly.Init() ;
    Matrix coeff = Tool::GetDoubleBlock(_constraints, 0, 0, _constraint_num,
        nonBasicNum) ;
    _glpk_poly.set_coefficients( std::move(coeff) ) ; 
    double currConstant ;
    for (int i = 0; i < _constraint_num; ++ i) {
      currConstant = Tool::GetDouble( (RNumber)_constraints.at(i,constraintsCols-1) ) ;
      _glpk_poly.SetConstant(i, currConstant) ;
    }
  }

  // objective in rational
  if ( ! removedVari.empty() ) {
    unsigned objRowSize = _objective.r.rows() ;
    unsigned objColSize = _objective.r.cols() ;
    RMatrix tmpObj(_objective.r) ;
    _objective.r.set_zero() ;
    for (unsigned i = 0; i < objRowSize; ++ i) {
      for (unsigned j = 0; j < remainedVari.size(); ++ j) {
        _objective.r.at(i,j) = tmpObj.at( i,remainedVari[j] ) ;
      }
      _objective.r.at(i,objColSize-1) = tmpObj.at(i,objColSize-1) ;
    }
  }

#ifdef DEBUGINFO_PLP
    log_mtx.lock() ;
    std::cout << "Objective: " << std::endl ;
    print(_objective.r) ;
    log_mtx.unlock() ;
#endif

  // objective in floating point
  double currF ;
  _objective.f.setZero() ;
  for (int i = 0; i < _parameter_num+1; ++ i) {
    for (int j = 0; j < _objective.r.cols(); ++ j) {
      currF = Tool::GetDouble( (RNumber)_objective.r.at(i, j) ) ;
      _objective.f(i, j) = currF ;
    }
  }

  int oriNonBasicNum = tmpConstraints.cols() - 1 ;
  if (newBdNum < oriNonBasicNum || ! removedVari.empty() ) {
    std::sort( bdIdx.begin(), bdIdx.end() ) ;
    _non_negative_idx.clear() ;
    std::vector<int> nonNegIdx(oriNonBasicNum, 0), reducedNegIdx ;
    for (int idx : bdIdx) {
      nonNegIdx[idx] = 1 ;
    }
    for (int i = 0; i < oriNonBasicNum; ++ i) {
      if ( std::find(removedVari.begin(), removedVari.end(), i) == removedVari.end() ) {
        reducedNegIdx.push_back( nonNegIdx[i] ) ;
      }
    }
    for (unsigned i = 0; i < reducedNegIdx.size(); ++ i) {
      if (reducedNegIdx[i] == 1) {
        _non_negative_idx.push_back(i) ;
      }
    } 
  }
  
#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "minimized constraints: " ;
  for (int i : consIdx) {
    std::cout << i << " " ;
  }
  std::cout << std::endl ;
  std::cout << "minimized bounds: " ;
  for (int i : _non_negative_idx) {
    std::cout << i << " " ;
  }
  std::cout << std::endl ;
  std::cout << "minimized constraints matrix: " ; print(_constraints) ;
  log_mtx.unlock() ;
#endif

  return true ;
}

Vector Plp::GetSimplexObj(const Vector& currPoint) {
  int nonBasicNum = _variable_num - _constraint_num ;
  Vector newObj(nonBasicNum) ;
  double sum, curr ;
  for (int j = 0; j < nonBasicNum; ++ j) {
    sum = 0.0 ;
    for (int i = 0; i < _parameter_num+1; ++ i) {
      curr = Tool::GetDouble( (RNumber)_objective.r.at(i, j) ) ;
      sum += curr * currPoint(i);
    }
    newObj(j) = sum ; 
  } 

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "glpk objective: " << std::endl << newObj << std::endl ;
  log_mtx.unlock() ;
#endif

  return newObj ;
}

RMatrix Plp::ReconstructCons(const std::vector<int>& basicIdx) const {
  return ReconstructCons(basicIdx, _constraints, true) ;
}

RMatrix Plp::ReconstructCons(const std::vector<int>& basicIdx,
    const RMatrix& cons, bool useDefSize) const {
  std::list<int> basisList ;
  for (unsigned i = 0; i < basicIdx.size(); ++ i) {
    basisList.push_back( basicIdx[i] ) ; 
  }
  int consNum ;
  int variNum = cons.cols() -1 ;
  if (useDefSize) {
    consNum = _constraint_num ;
  }
  else {
    consNum = cons.rows() ;
  }
  RNumber den, curr, ratio ;
  RMatrix coeff(cons) ;
  std::list<int>::iterator it ;
  int currIdx, rowIdx ;

  std::vector<int> rowOrder ;

  for (int i = 0; i < consNum; ++ i) {
    currIdx = -1 ;
    it = basisList.begin() ;
    while (it != basisList.end() ) {
      if ( ! coeff.at(i, *it).is_zero() ) {
        currIdx = *it ;
        basisList.erase(it) ;
        break ;
      }
      ++ it ; 
    }
    if (currIdx == -1) {
      //std::cerr << "ReconstructCons: cannot find basis." << std::endl ;
      //std::terminate() ;
      return RMatrix(0, 0) ;
    }

    auto begin = basicIdx.begin() ;
    auto current = std::find(basicIdx.begin(), basicIdx.end(), *it) ;

    rowIdx = std::distance(begin, current) ; 
    rowOrder.push_back(rowIdx) ;
    den = coeff.at(i, currIdx) ;
    for (int j = 0; j < variNum+1; ++ j) {
      curr = coeff.at(i, j) ;
      coeff.at(i, j) = curr/den ; 
    }
    for (int k = 0; k < (int)basicIdx.size(); ++ k) {
      if (k == i) continue ;
      ratio = coeff.at(k, currIdx) ;
      for (int j = 0; j < variNum+1; ++ j) {
        curr = coeff.at(k, j) ;
        coeff.at(k, j) = curr - ratio * coeff.at(i, j) ;
      } 
    }
  }

  RMatrix orderedCoeff( consNum, coeff.cols() ) ;
  for (int i = 0; i < consNum; ++ i) {
    for (unsigned j = 0; j < coeff.cols(); ++ j) {
      orderedCoeff.at(rowOrder[i],j) = coeff.at(i,j) ;
    }
  }

  return orderedCoeff ;
}

bool Plp::VerifyGlpkResult(const std::vector<int>& basis) {
  RMatrix coeff( ReconstructCons(basis) ) ; 
  RNumber zero ;
  int constIdx = coeff.cols() - 1 ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "Verify result of GLPK: constraints: " << std::endl ; 
  print(coeff) ;
  log_mtx.unlock() ;
#endif

  for (int i = 0; i < _constraint_num; ++ i) {
    if (coeff.at(i, constIdx) < zero) {
      return false ;
    }
  } 

  return true ;
}

/*******************************************************************************
 * The constraints are stored as row vector: ax + b >= 0
*******************************************************************************/
bool Plp::ConsRedundantR(const RMatrix& matrix, int idx) {
  int consNum = matrix.rows() ;
  int colNum = matrix.cols() ;
  RMatrix farkasM( matrix.transpose() ) ;
  for (int j = 0; j < colNum; ++ j) {
    farkasM.at(j,idx) = matrix.at(consNum-1,j) ;
    farkasM.at(j,consNum-1) = matrix.at(idx,j) ;
  }

  ReducedMatrix echelonM = Tool::GaussianElimination(farkasM) ;
  // if the problem is infeasible, e.g. 0 = -2
  if ( ! echelonM.succeed ) {
    return false ; 
  }

  std::vector<int> rowOrder = echelonM.rowIdx ;
  int rowSize = rowOrder.size() ;
  int colSize = echelonM.colIdx.size() ;
  int initColSize = echelonM.matrix.cols() ;
  std::vector<int> basisOrder = echelonM.colIdx ;
  RMatrix reducedM(rowSize, colSize+1) ;
  for (int i = 0; i < rowSize; ++ i) {
    for (int j = 0; j < colSize; ++ j) {
      reducedM.at(i, j) = echelonM.matrix.at( rowOrder[i], basisOrder[j] ) ;
    }
    reducedM.at(i, colSize) = echelonM.matrix.at(rowOrder[i], initColSize-1) ;
  }

  Simplex solver(reducedM, rowSize, colSize) ;
  return ! solver.Solve(true).empty() ;
}

std::vector<int> Plp::MinimizePolyR(const RMatrix& matrix) {
  int consNum = matrix.rows() ;
  std::vector<int> activeIdx ;
  for (int k = 0; k < consNum; ++ k) {
    if( ! ConsRedundantR(matrix, k) ) {
      activeIdx.push_back(k) ;
    }
  } 

  return activeIdx ;
}

/*******************************************************************************
 * The constraints are stored as row vector: ax + b >= 0
*******************************************************************************/
bool Plp::IsFlatRegionR(const RMatrix& matrix) {
  int rowNum = matrix.rows() ;
  int colNum = matrix.cols() ;
  int newColNum = 2 * colNum + rowNum - 1 ;
  // the parameters are not non-negative,
  // construct the constraints in the form of -a(x_1 - x_2) + r = b-shift
  RMatrix shiftM(rowNum, newColNum) ;
  RNumber one ;
  one.set_one() ;
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum-1; ++ j) {
      shiftM.at(i,2*j) = matrix.at(i,j) ;
      shiftM.at(i,2*j+1) = - matrix.at(i,j) ;
    }
    shiftM.at(i,i+2*colNum-2) = one ;
    shiftM.at(i,newColNum-1) = matrix.at(i,colNum-1) - one ; 
  }
  
  ReducedMatrix echelonM = Tool::GaussianElimination(shiftM) ;
  // if the problem is infeasible, e.g. 0 = -2
  if ( ! echelonM.succeed) {
    return true ;
  }

  std::vector<int> rowOrder = echelonM.rowIdx ;
  std::vector<int> colOrder = echelonM.colIdx ;
  int rowSize = rowOrder.size() ;
  int colSize = echelonM.colIdx.size() ;
  int initColSize = echelonM.matrix.cols() ;
  std::vector<int> basisOrder = echelonM.colIdx ;
  RMatrix reducedM(rowSize, colSize+1) ;
  for (int i = 0; i < rowSize; ++ i) {
    for (int j = 0; j < colSize; ++ j) {
      reducedM.at(i, j) = echelonM.matrix.at( rowOrder[i], basisOrder[j] ) ;
    }
    reducedM.at(i, colSize) = echelonM.matrix.at(rowOrder[i], initColSize-1) ;
  }

  Simplex solver(reducedM, rowSize, colSize) ;
  std::vector<int> basicIdx = solver.Solve(true) ;
  if( basicIdx.empty() ) {
    return true ;
  }
  else {
    /*
    RMatrix feasibleRes( solver.get_feasible_vertex() ) ;
    RNumber zero ;
    std::vector<RNumber> variVal(newColNum-1, zero) ;
    assert( newColNum-1 == colOrder.size() ) ;
    for (unsigned i = 0; i < basicIdx.size(); ++ i) {
      variVal[ colOrder[ basicIdx[i] ] ] = feasibleRes.at(0, i) ;
    }
    RMatrix point(1,colNum-1) ;
    for (int i = 0; i < colNum-1; ++ i) {
      point.at(0,i) = variVal[2*i] - variVal[2*i+1] ;
    }
    */
    return false ;
  }
}

#ifdef VERIFY_PLP
bool Plp::VerifyCentralPoint(const RMatrix& point, const Polyhedron& poly) {
  int consNum = poly.get_constraint_num() ;
  int variNum = poly.get_variable_num() ;
  RMatrix rationalPoly(consNum, variNum) ;
  RMatrix rationalConstant(consNum, 1) ;
  long int val ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      val = poly.GetCoef(i, j) ;
      rationalPoly.at(i, j) = val ; 
    }
    val = poly.GetConstant(i) ;
    rationalConstant.at(i, 0) = val ;
  }
  RMatrix result(rationalPoly * point - rationalConstant) ;
  RNumber zero ;
  for (int i = 0; i < consNum; ++ i) {
    if (result.at(i, 0) >= zero) {
      return false ;
    }
  }
  return true ;
}
#endif
