/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
*******************************************************************************/

#include <iostream>
#include <utility>
#include <random>
#include <cmath>
#include <fstream>
#include "pplp_raytracing.h"
#include "pplp_apronInterface.h"
#include "pplp_glpkInterface.h"

using namespace PPLP ;

Polyhedron::Polyhedron() 
  : _constraint_num(0), _eq_constraint_num(0), _variable_num(0),
  _redundant_num(-1), _zero_num(-1),
  _generator_num(-1), _id(-1), _is_minimized(false), _status(PolyStatus::top) 
{}

void Polyhedron::SetSize(int consNum, int variNum, int eqConsNum) {
  _constraint_num = consNum ;
  _variable_num = variNum ; 
  _eq_constraint_num = eqConsNum ;
  _coefficients.resize(consNum, variNum) ;
  _constants.resize(consNum) ;
  _eq_coefficients.resize(eqConsNum, variNum) ;
  _eq_constants.resize(eqConsNum) ;
  _coefficients.setZero() ;
  _constants.setZero() ;
  _eq_coefficients.setZero() ;
  _eq_constants.setZero() ;
  _active_table.resize(consNum, ConstraintState::redundant) ;
  _witness_point.resize(consNum) ;
  if (consNum == 0 && eqConsNum == 0) {
    _status = PolyStatus::top ;
  }
  else {
    _status = PolyStatus::normal ;
  }
}

Polyhedron::Polyhedron(int consNum, int variNum, int eqConsNum) 
  : _constraint_num(consNum), _eq_constraint_num(eqConsNum), _variable_num(variNum),
  _redundant_num(-1), _zero_num(-1), _generator_num(-1),
  _id(-1), _is_minimized(false), _status(PolyStatus::normal) {
  _coefficients.resize(consNum, variNum) ;
  _constants.resize(consNum) ;
  _eq_coefficients.resize(eqConsNum, variNum) ;
  _eq_constants.resize(eqConsNum) ;
  _coefficients.setZero() ;
  _constants.setZero() ;
  _eq_coefficients.setZero() ;
  _eq_constants.setZero() ;
  _active_table.resize(consNum, ConstraintState::redundant) ;
  _witness_point.resize(consNum) ;
  _status = PolyStatus::normal ;
}

Polyhedron::Polyhedron(const RMatrix& matrix, const std::vector<int>& idx, const RMatrix& eqMatrix, const std::vector<int>& eqIdx) 
  : _redundant_num(-1), _zero_num(-1), _generator_num(-1), _id(-1),
  _is_minimized(false), _status(PolyStatus::normal) {
  std::vector<int> idxCopy = idx ;
  if ( idxCopy.empty() ) {
    for (int i = 0; i < matrix.rows(); ++ i) {
      idxCopy.push_back(i) ;
    }
  }
  std::vector<int> eqIdxCopy = eqIdx ;
  if ( eqIdxCopy.empty() ) {
    for (int i = 0; i < eqMatrix.rows(); ++ i) {
      eqIdxCopy.push_back(i) ;
    }
  }

  int consNum = idxCopy.size() ;
  int variNum = matrix.cols() - 1 ; 
  int eqConsNum = eqIdxCopy.size() ;
  SetSize(consNum, variNum, eqConsNum) ;
  int currIdx ;
  for (int i = 0; i < consNum; ++ i) {
    currIdx = idxCopy[i] ;
    for (int j = 0; j < variNum; ++ j) {
      _coefficients(i,j) = Tool::GetDouble( (RNumber)matrix.at(currIdx,j) ) ;
    }
    _constants(i) = Tool::GetDouble( (RNumber)matrix.at(currIdx,variNum) ) ;
  }

  // equalities
  for (int i = 0; i < eqConsNum; ++ i) {
    currIdx = eqIdxCopy[i] ;
    for (int j = 0; j < variNum; ++ j) {
      _eq_coefficients(i,j) = Tool::GetDouble( (RNumber)eqMatrix.at(currIdx,j) ) ;
    }
    _eq_constants(i) = Tool::GetDouble( (RNumber)eqMatrix.at(currIdx,variNum) ) ;
  } 
}

/*******************************************************************************
 * Set the polyhedron to initial state. You can reset the constraints, but the
 * original ones will not be deleted.
*******************************************************************************/
void Polyhedron::Init() {
  _is_minimized = false ;
  if (_constraint_num == 0 && _eq_constraint_num == 0) {
    _status = PolyStatus::top ;
  }
  else {
    _status = PolyStatus::normal ;
  }
  _coefficients.setZero() ;
  _constants.setZero() ;
  _eq_coefficients.setZero() ;
  _eq_constants.setZero() ;
  _active_table.clear() ;
  _active_table.resize(_constraint_num, ConstraintState::redundant) ;  
  _internal_point.Clear() ;
  _witness_point.clear() ;
  _witness_point.resize(_constraint_num) ;
  _redundant_num = -1 ;
  _zero_num = -1 ;
  _generator_num = -1 ;
  _id = -1 ;
}

/*******************************************************************************
 * Substitute equalities into inequalities 
 * If the coefficients are integers (the input polyhedra always have integers as
 * the coefficients), we use rational numbers in this function.
 * Otherwise we use floating point numbers and a threshold, but this cannot
 * guarantee soundness.
 * In projection the constraints of the regions are all inequalities, so there
 * will be no impact. 
*******************************************************************************/
ReducedMatrixFloat Polyhedron::SubsEqConsF() {
  if (_eq_constraint_num == 0) {
    return ReducedMatrixFloat() ; 
  }

  Matrix redM ;
  int rowNum, variNum, rIdx ;
  Matrix m(_eq_constraint_num, _variable_num+1) ;
  m.block(0, 0, _eq_constraint_num, _variable_num) = _eq_coefficients ;
  m.col(_variable_num) = _eq_constants.transpose() ;
  ReducedMatrixFloat redMatrixF = Tool::GaussianEliminationFloat(m) ;
  if ( ! redMatrixF.succeed) {
    return redMatrixF ;
  }

  rowNum = redMatrixF.rowIdx.size() ;
  variNum = redMatrixF.colIdx.size() ;
  redM.resize(rowNum, variNum+1) ;
  _eq_coefficients.resize(rowNum, _variable_num) ;
  _eq_constants.resize(rowNum) ;
  for (int i = 0; i < rowNum; ++ i) {
    rIdx = redMatrixF.rowIdx[i] ;
    _eq_coefficients.row(i) = redMatrixF.matrix.row(rIdx).head(_variable_num) ;
    _eq_constants(i) = redMatrixF.matrix(rIdx, _variable_num) ;
  }
  _eq_constraint_num = rowNum ;
  return redMatrixF ;
}

ReducedMatrix Polyhedron::SubsEqConsR() {
  if (_eq_constraint_num == 0) {
    return ReducedMatrix(0,0) ; 
  }

  Matrix redM ;
  int rowNum, variNum, rIdx, currVal ;
  RMatrix rmat(_eq_constraint_num, _variable_num+1) ;
  for (int i = 0; i < _eq_constraint_num; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      currVal = _eq_coefficients(i,j) ;
      rmat.at(i,j) = currVal ;
    }
    currVal = _eq_constants(i) ;
    rmat.at(i,_variable_num) = currVal ;
  }

  ReducedMatrix redMatrix = Tool::GaussianElimination(rmat) ; 
  if ( ! redMatrix.succeed) {
    return redMatrix ;
  }

  rowNum = redMatrix.rowIdx.size() ;
  variNum = redMatrix.colIdx.size() ;
  redM.resize(rowNum, variNum+1) ;
  _eq_coefficients.resize(rowNum, _variable_num) ;
  _eq_constants.resize(rowNum) ;
  for (int i = 0; i < rowNum; ++ i) {
    rIdx = redMatrix.rowIdx[i] ;
    for (int j = 0; j < variNum; ++ j) {
      _eq_coefficients(i, j) = Tool::GetDouble( (RNumber)redMatrix.matrix.at(rIdx, j) ) ;

    }
    _eq_constants(i) = Tool::GetDouble( (RNumber)redMatrix.matrix.at(rIdx, _variable_num) ) ;

  }
  _eq_constraint_num = rowNum ;

  return redMatrix ;
}


void Polyhedron::SubsIneqConsF(const ReducedMatrixFloat& redM) {
  // inequalities
  Matrix ineqFmat(_constraint_num, _variable_num+1) ;
  ineqFmat.block(0, 0, _constraint_num, _variable_num) = _coefficients ;
  ineqFmat.col(_variable_num) = _constants.transpose() ;

  double curr, ratio ;
  int basicIdx, eqRowIdx ;
  for (int i = 0; i < _eq_constraint_num; ++ i) {
    eqRowIdx = redM.rowIdx[i] ;
    basicIdx = redM.basicIdx[i] ;
    for (int k = 0; k < _constraint_num; ++ k) {
      ratio = ineqFmat(k, basicIdx) ;
      for (int j = 0; j < _variable_num+1; ++ j) {
        curr = ineqFmat(k,j) - ratio * redM.matrix(eqRowIdx, j) ;
        ineqFmat(k,j) = (curr < 1e-6 && curr > -1e-6) ? 0 : curr ;  
      }
    }
  }

  std::vector<int> nonDupIdx = Tool::GetNonDupIdx(ineqFmat, true) ;
  _coefficients.resize(nonDupIdx.size(), _variable_num) ;
  _constants.resize( nonDupIdx.size() ) ;
  for (unsigned i = 0; i < nonDupIdx.size(); ++ i) {
    _coefficients.row(i) = ineqFmat.row(i).head(_variable_num) ;
    _constants(i) = ineqFmat(i, _variable_num) ;
  }
  _constraint_num = nonDupIdx.size() ;

}

void Polyhedron::SubsIneqConsR(const ReducedMatrix& redM) {
  int currVal ;
  // inequalities
  RMatrix ineqRmat(_constraint_num, _variable_num+1) ;
  for (int i = 0; i < _constraint_num; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      currVal = _coefficients(i,j) ;
      ineqRmat.at(i,j) = currVal ;
    }
    currVal = _constants(i) ;
    ineqRmat.at(i,_variable_num) = currVal ;
  }

  RNumber currR ;
  int basicIdx, eqRowIdx ;
  for (int i = 0; i < _eq_constraint_num; ++ i) {
    eqRowIdx = redM.rowIdx[i] ;
    basicIdx = redM.basicIdx[i] ;
    for (int k = 0; k < _constraint_num; ++ k) {
      currR = ineqRmat.at(k, basicIdx) ;
      for (int j = 0; j < _variable_num+1; ++ j) {
        ineqRmat.at(k,j) = ineqRmat.at(k,j) - currR * redM.matrix.at(eqRowIdx, j) ;
      }
    }
  }

  std::vector<int> nonDupIdx = Tool::GetNonDupIdx(ineqRmat) ;
  _coefficients.resize(nonDupIdx.size(), _variable_num) ;
  _constants.resize( nonDupIdx.size() ) ;
  int rIdx ;
  for (unsigned i = 0; i < nonDupIdx.size(); ++ i) {
    rIdx = nonDupIdx[i] ;
    for (int j = 0; j < _variable_num; ++ j) {
      _coefficients(i,j) = Tool::GetDouble( (RNumber)ineqRmat.at(rIdx, j) ) ;
    }
    _constants(i) = Tool::GetDouble( (RNumber)ineqRmat.at(rIdx, _variable_num) ) ;
  }
  _constraint_num = nonDupIdx.size() ;
}


/*******************************************************************************
 * Calls the raytracing method to eliminate the redundant constraints, and sets
 * the corresponding value in _active_table as false
 * @return true if minimize successful
*******************************************************************************/
bool Polyhedron::Minimize(bool getWitness, bool cone, bool useFloat) {

#ifdef DEBUGINFO_RAYTRACING
  log_mtx_raytracing.lock() ; 
  std::cout << "Raytracing minimization starts." << std::endl ;
  log_mtx_raytracing.unlock() ; 
#endif

  if ( IsTop() ) {
    std::cerr << "Minimization failed. There is no constraints." << std::endl ;
    std::terminate() ;
  }
  if (_is_minimized) {
    return true ;
  }

  // start to minimize
  // substitute equalities into inequalities
  if (_eq_constraint_num != 0) {
    if (useFloat) {
      ReducedMatrixFloat redMatrixF = SubsEqConsF() ;
      if ( ! redMatrixF.succeed) {
        SetBottom() ;
        return false ;
      }
      if (_constraint_num != 0) {
        SubsIneqConsF(redMatrixF) ;
      }
    }
    else {
      ReducedMatrix redMatrixR = SubsEqConsR() ;
      if ( ! redMatrixR.succeed) {
        SetBottom() ;
        return false ;
      }
      if (_constraint_num != 0) {
        SubsIneqConsR(redMatrixR) ;
      }
    }
  }

  if (_constraint_num == 0) {
    return true ;
  }
  // raytracing
  if ( _internal_point.IsEmpty() ) {
    set_internal_point( GlpkInterface::GetCentralPoint(*this) ) ;
    // for unbounded polyhedron
    if ( _internal_point.IsEmpty() ) {
      set_internal_point( GlpkInterface::GetSatPoint(*this) ) ;
    }
    if(_internal_point.IsEmpty()) {
      std::cerr << "Minimization failed. Cannot find a point inside the polyhedra" << std::endl ;
      return false ;
    }
  }

#ifdef DEBUGINFO_RAYTRACING
  log_mtx_raytracing.lock() ; 
  std::cout << "Found interior point: " << _internal_point << std::endl ;
  log_mtx_raytracing.unlock() ; 
#endif

  Raytracing raytrace(*this, _internal_point, getWitness, cone) ;
  raytrace.RayHitting() ; 
  raytrace.Determine() ;
  _is_minimized = true ;
  _redundant_num = GetInactiveIdx().size() ;

  return true ;
}

/*******************************************************************************
 * Minimization without raytracing.
 * TODO this function should be abandoned or remove the threshod in Sat()
*******************************************************************************/
void Polyhedron::MinimizeSimple() {
  std::vector<bool> res(_constraint_num) ;
  // #pragma omp parallel for
  for (int i = 0; i < _constraint_num; ++ i) {
    res[i] = GlpkInterface::Sat(*this, i) ;
    if (res[i] == true) {
      Activate(i) ;       
    }
  }  
  _is_minimized = true ;
}

/*******************************************************************************
 * Gets the minimized polyhedron. If the polyhedron is not minimized, minimize
 * it first. Copy all the reasonal data of the original polyhedron.
 * @return the minimized polyhedron
*******************************************************************************/
Polyhedron Polyhedron::GetMinimizedPoly(bool getWitness, bool cone) {
  if (! _is_minimized) {
    Minimize(getWitness, cone) ;
  }
  std::vector<int> indexVec = GetActiveIdx() ; 
  Polyhedron miniPoly = GetSubPoly(indexVec) ;
  miniPoly.set_internal_point(_internal_point) ;
  int consNum = miniPoly.get_constraint_num() ;
  for (int i = 0; i < consNum; ++ i) {
    miniPoly._active_table[i] = ConstraintState::irredundant ;
    miniPoly._witness_point[i] = _witness_point[ indexVec[i] ] ;
  }
  miniPoly._is_minimized = true ;
  return miniPoly ;
}

/*******************************************************************************
 * Creates a polyhedron which contains a subset of constraints of the current
 * polyhedron. Just copy the constraints.
 * @para index the index of the constraints in the subset
 * @return the new polyhedron contains the subset of constraints 
*******************************************************************************/
Polyhedron Polyhedron::GetSubPoly(const std::vector<int>& indexVec, const std::vector<int>& eqIndexVec) {
  std::vector<int> idxCopy, eqIdxCopy ;
  if ( indexVec.empty() ) {
    for (int i = 0; i < _constraint_num; ++ i) {
      idxCopy.push_back(i) ;
    } 
  }
  else {
    idxCopy = indexVec ;
  }

  if ( eqIndexVec.empty() ) {
    for (int i = 0; i < _eq_constraint_num; ++ i) {
      eqIdxCopy.push_back(i) ;
    }
  }
  else {
    eqIdxCopy = eqIndexVec ;
  }

  int consNum = idxCopy.size() ;
  int variNum = get_variable_num() ;
  int eqConsNum = eqIdxCopy.size() ;
  Polyhedron newPoly(consNum, variNum, eqConsNum) ;
  int currIdx ;
  for (int i = 0; i < consNum; ++ i) {
    currIdx = idxCopy[i] ;
    for (int j = 0; j < variNum; ++ j) {
      newPoly.SetCoef( i, j, GetCoef(currIdx, j) ) ;
    }
    newPoly.SetConstant( i, GetConstant(currIdx) ) ;
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    currIdx = eqIdxCopy[i] ;
    for (int j = 0; j < variNum; ++ j) {
      newPoly.SetEqCoef( i, j, GetEqCoef(currIdx, j) ) ;
    }
    newPoly.SetEqConstant( i, GetEqConstant(currIdx) ) ;
  }
  newPoly._zero_num = _zero_num ;
  return newPoly ;
}

/*******************************************************************************
 * Compare two constraints in lexico order. 
 * cons1 and cons2 containts the coefficients of variables and also the constant
 * @return CmpOperator::greater if cons1 > cons2 
*******************************************************************************/
CmpOperator Polyhedron::CmpConstraint(const Vector& cons1, const Vector& cons2) {
  CmpOperator res = CmpOperator::equal ;
  for (int j = 0; j < cons1.size(); ++ j) {
    if ( cons1(j) > cons2(j) ) {
      res = CmpOperator::greater ;
      break ;
    }
    else if ( cons1(j) < cons2(j) ) {
      res = CmpOperator::less ;
      break ;
    }
  }
  return res ;
}

/*******************************************************************************
 * Selection sorting 
*******************************************************************************/
void Polyhedron::SortMatrix(Matrix& m) {
  int consNum = m.rows() ;
  int minIdx ;
  Vector tmpVec ;
  for (int i = 0; i < consNum; ++ i) {
    minIdx = i ;
    for (int k = i+1; k < consNum; ++ k) { 
      if (CmpConstraint( m.row(minIdx), m.row(k) ) == CmpOperator::greater) {
        minIdx = k ;
      } 
    }
    if (minIdx != i) {
      tmpVec = m.row(i) ;
      m.row(i) = m.row(minIdx) ;
      m.row(minIdx) = tmpVec ;
    }
  } 
}

/*******************************************************************************
 * Checks equality by sorting in lexico order. 
 * @para poly the polyhedron to compare
 * @return true if the two polyhedra are equal
*******************************************************************************/
bool Polyhedron::IsEqualTo(const Polyhedron& poly) {
  if ( get_coefficients().size() != poly.get_coefficients().size() ) {
    return false ;
  }
  if ( get_eq_coefficients().size() != poly.get_eq_coefficients().size() ) {
    return false ;
  }

  int consNum = poly.get_constraint_num() ;
  int variNum = poly.get_variable_num() ;
  Matrix matrix1(consNum, variNum+1), matrix2(consNum, variNum+1) ;
  matrix1.block(0, 0, consNum, variNum) = get_coefficients() ;
  matrix1.col(variNum) = get_constants().transpose() ;
  matrix2.block(0, 0, consNum, variNum) = poly.get_coefficients() ;
  matrix2.col(variNum) = poly.get_constants().transpose() ;
  SortMatrix(matrix1) ; 
  SortMatrix(matrix2) ;
  bool res1 = matrix1 == matrix2 ;

  int eqConsNum = poly.get_eq_constraint_num() ;
  Matrix eqMatrix1(eqConsNum, variNum+1), eqMatrix2(eqConsNum, variNum+1) ;
  eqMatrix1.block(0, 0, eqConsNum, variNum) = get_eq_coefficients() ;
  eqMatrix1.col(variNum) = get_eq_constants().transpose() ;
  eqMatrix2.block(0, 0, eqConsNum, variNum) = poly.get_eq_coefficients() ;
  eqMatrix2.col(variNum) = poly.get_eq_constants().transpose() ;
  SortMatrix(eqMatrix1) ; 
  SortMatrix(eqMatrix2) ;
  bool res2 = eqMatrix1 == eqMatrix2 ;
  
  return res1 && res2 ;
}

/*******************************************************************************
 * Checks inclusion with raytracing. 
 * For example we check if p1 is included in P2, we need three steps:
 * 1) Check where is the inner point x of P1. If x is not inside P2, then P1
 * is NOT included in P2; 
 * 2) Check each constraint in P2. If any constraint in P2 is irredundant wrt P1
 * then P1 is not included in P2;
 * 3) If no constraint in P2 is irredundant wrt P1 then P1 is included in P2.
 * Note that P1 and P2 should be minimized.
 * @para poly the polyhedron which may be included in the current one, i.e. here
 * we check if poly is included in *this.
 * @return true if poly is included in the current polyhedron 
*******************************************************************************/
bool Polyhedron::Include(Polyhedron& poly) {
  if ( get_variable_num() != poly.get_variable_num() ) {
    std::cerr << "Caonnt check inclusion. Polyhedra are not in the same environment." 
      << std::endl ;
    std::terminate() ;
  }
  const int consNumP1 = poly.get_constraint_num() ;
  const int consNumP2 = get_constraint_num() ;
  if ( poly.get_internal_point().IsEmpty() ) {
    poly.set_internal_point( GlpkInterface::GetCentralPoint(poly) ) ;
  }

  // check if x is inside P2
  for (int i = 0; i < consNumP2; ++ i) {
    if (Satisfy(poly.get_internal_point(), i) == false) {
      return false ;
    }
  } 
  Polyhedron newPoly = Combine(poly) ;
  newPoly.set_internal_point( poly.get_internal_point() ) ;
  if(newPoly._internal_point.IsEmpty()) {
    std::cerr << "Cannot check inclusion. Cannot find a point inside the polyhedra" 
      << std::endl ;
    std::terminate() ;
  }
  Raytracing raytrace(newPoly, newPoly._internal_point, false, false) ;
  // check constraints of P2  
  raytrace.SetInclusion(consNumP1) ; 
  raytrace.RayHitting() ;    
  if (raytrace.HasInclusion() == false) {
    return false ;
  }
  raytrace.Determine() ;
  if (raytrace.HasInclusion() == false) {
    return false ;
  }
  
  return true ;
}

/*******************************************************************************
 * Checks inclusion with Farkas's lemma. 
 * @para poly the polyhedron which may be included in the current one, i.e. here
 * we check if poly is included in *this.
 * @return true if poly is included in the current polyhedron 
*******************************************************************************/
bool Polyhedron::IncludeStandard(Polyhedron& poly) {
  int consNum = poly.get_constraint_num() ;
  int variNum = poly.get_variable_num() ;
  for (int i = 0; i < _constraint_num; ++ i) {
    // constraints in matrix are Cx + b >= 0
    // aij, bi are coeff of poly, pij,qi are coeff of *this, the LP is:
    // a11 ... an1 0             pj1
    // ... ... ... .   lambda =  ...
    // a1n ... ann 0             pjn
    // b1  ...  bn 1             qj
    Polyhedron simplexPoly(0, consNum+1, variNum+1) ;
    simplexPoly._eq_coefficients.block(0, 0, variNum, consNum) = 
        - poly.get_coefficients().transpose() ; 
    simplexPoly._eq_coefficients.row(variNum).head(consNum) = poly.get_constants() ;
    simplexPoly._eq_coefficients.col(consNum).head(variNum).setZero() ;
    simplexPoly._eq_coefficients(variNum, consNum) = 1.0 ;
    simplexPoly._eq_constants.head(variNum) = - get_coefficients().row(i).transpose() ;
    simplexPoly._eq_constants(variNum) = GetConstant(i) ;

    /*
    for (int p = 0; p < variNum+1; ++ p) {
      simplexPoly.SetOperator(p, ConstraintOperator::equal) ;
    }
    */

    GlpkInterface glpkInter ;
    bool simplexRes = glpkInter.Simplex(simplexPoly, Vector(), true, true, false) ;
    if (simplexRes == false) {
      return false ;
    }
  }
  return true ;
}

/*******************************************************************************
 * Put constraints in two polyhedra together in a new polyhedron. The polyhedron
 * at the parameter is put at the beginning.
 * @para poly the polyhedron to be put at the beginning
 * @return the new polyhedron contains all constraints 
*******************************************************************************/
Polyhedron Polyhedron::Combine(const Polyhedron& poly) {
  int consNum1 = poly.get_constraint_num() ;
  int consNum2 = get_constraint_num() ;
  int variNum = get_variable_num() ;
  Polyhedron newPoly(consNum1+consNum2, variNum) ;
  newPoly._coefficients.block(0, 0, consNum1, variNum) = poly.get_coefficients() ;
  newPoly._coefficients.block(consNum1, 0, consNum2, variNum) = get_coefficients() ;
  newPoly._constants.head(consNum1) = poly.get_constants() ;
  newPoly._constants.tail(consNum2) = get_constants() ;
  return newPoly ;
}

/*******************************************************************************
 * Checks if a point satisfies the constraint.
 * @para point the point to be satisfied 
 * @para index the index of the constraint to be checked
 * @para strict true when the point should not be on the boundary
 * @return true if the point satisfies the constraint
*******************************************************************************/
bool Polyhedron::Satisfy(const Point& point, int index, bool strict) const {
  if ( point.IsEmpty() ) {
    std::cerr << "Satisfy() failed. Point is empty." << std::endl ;
    std::terminate() ;
  }
  double res = point.get_coordinates() * get_coefficients().row(index).transpose()  ;
  double threshold = Tool::GetDotProductThreshold( point.get_coordinates(),
      get_coefficients().row(index) ) ;

#ifdef DEBUGINFO_RAYTRACING
  log_mtx_raytracing.lock() ; 
  std::cout << "test: dot product: " << res << " const: " << GetConstant(index) << std::endl ;
  log_mtx_raytracing.unlock() ; 
#endif

  if (strict) {
    if( Double::IsLessThan(res, GetConstant(index), threshold) ) {
      return true ;
    }
  } 
  else {
    if( Double::IsLessEq(res, GetConstant(index), threshold) ) {
      return true ;
    }
  }
  return false ;
}

/*******************************************************************************
 * Checks if a point is inside a polyhedron.
 * @para point the point to be satisfied 
 * @para strict true when the point should not be on the boundary
 * @return true if the point is inside the polyhedron
*******************************************************************************/
bool Polyhedron::Satisfy(const Point& point, bool strict) const {
  if ( point.IsEmpty() ) {
    std::cerr << "Satisfy() failed. Point is empty." << std::endl ;
    std::terminate() ;
  }
  bool sat = true ;
  for (int i = 0; i < get_constraint_num(); ++ i) {
    sat = Satisfy(point, i, strict) ;
    if (sat == false) break ;
  }
  return sat ;
}

/*******************************************************************************
 * Checks if a point satisfies the constraint.
 * @para point the point to be satisfied 
 * @para index the index of the constraint to be checked
 * @para strict true when the point should not be on the boundary
 * @return the index of the constraint if the point is on the boundary,
 * otherwise returns inside or outside 
*******************************************************************************/
int Polyhedron::PointOnBoundary(const Point& point) const {
  if ( point.IsEmpty() ) {
    std::cerr << "Satisfy() failed. Point is empty." << std::endl ;
    std::terminate() ;
  }
  int sat = PointPos::inside ;
  double res, threshold ;
  for (int i = 0; i < get_constraint_num(); ++ i) {
    res = point.get_coordinates() * get_coefficients().row(i).transpose()  ;
    threshold = Tool::GetDotProductThreshold( point.get_coordinates(),
        get_coefficients().row(i) ) ;
    if ( Double::AreEqual(res, GetConstant(i), threshold) ) {
      sat = i ;
      // don't break
    }
    if ( Double::IsLessThan(GetConstant(i), res, threshold) ) {
      return PointPos::outside ;
    }
  }
  return sat ;
}

std::vector<int> Polyhedron::GetActiveIdx() const {
  if (! _is_minimized) {
    std::cerr << "Cannot get the active index. The polyhedron is not minimized." 
        << std::endl ;
    std::terminate() ;
  }
  std::vector<int> idxVec ;
  for (int i = 0; i < (int)_active_table.size(); ++ i) {
    if (_active_table[i] == ConstraintState::irredundant) {
      idxVec.push_back(i) ; 
    } 
  }
  return idxVec ;
}

std::vector<int> Polyhedron::GetInactiveIdx() const {
  if (! _is_minimized) {
    std::cerr << "Cannot get the inactive index. The polyhedron is not minimized." 
        << std::endl ;
    std::terminate() ;
  }
  std::vector<int> idxVec ;
  for (int i = 0; i < (int)_active_table.size(); ++ i) {
    if (_active_table[i] != ConstraintState::irredundant) {
      idxVec.push_back(i) ; 
    } 
  }
  return idxVec ;
}

/*******************************************************************************
 * _active_table[i] is 1 if it is irredundant, 0 if it is redundant.
 * @return true if the constrinat is irredundant
 * Note: this functions may be not thread safe
*******************************************************************************/
bool Polyhedron::IsActive(int consIdx) const {
  bool ret ;
  if (_active_table[consIdx] == ConstraintState::irredundant) {
    ret = true ;
  }
  else
  {
    ret = false ;
  }
  return ret ;
}

/*
std::string Polyhedron::GetOperatorStr(int consIdx) const {
  ConstraintOperator op = _operators[consIdx] ;
  if (op == ConstraintOperator::less) return "<" ;
  else if (op == ConstraintOperator::lesseq) return "<=" ;
  else if (op == ConstraintOperator::equal) return "=" ;
  else return "" ;
}
*/

/*******************************************************************************
 * @return true if the polyhedron is open
*******************************************************************************/
bool Polyhedron::IsOpen() const {
  Vector obj(_variable_num) ;
  obj.setOnes() ;
  GlpkInterface glpkInter1, glpkInter2 ;
  bool res1 = glpkInter1.Simplex(*this, -obj, false, false, false) ;
  bool res2 = glpkInter2.Simplex(*this, obj, false, false, false) ;
  return ! (res1 && res2) ;
}

/*******************************************************************************
 * @return true if the polyhedron is minimized
*******************************************************************************/
bool Polyhedron::IsMinimized() const {
  return _is_minimized ;
}

/*******************************************************************************
 * Prints the index of the avtivate constraints
*******************************************************************************/
void Polyhedron::PrintActiveIdx() const {
  std::cout<< "The index of the activate constraints are: " << std::endl ;
  for (int i = 0; i < (int)_active_table.size(); ++ i) {
    if (_active_table[i] == ConstraintState::irredundant) {
      std::cout << i << " " ;
    } 
  }
  std::cout << std::endl ;
}

/*******************************************************************************
 * Prints the constraints to the terminal
*******************************************************************************/
void Polyhedron::Print() const {
  std::cout << _constraint_num + _eq_constraint_num
      << " constraints in " << _variable_num  
      << " dimension (in form [coefficients <op> constant]):" << std::endl ;
  for (int i = 0; i < _eq_constraint_num; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      std::cout << _eq_coefficients(i, j) << " " ;
    }
    std::cout << " = " << _eq_constants[i] << std::endl ;
  }
  for (int i = 0; i < _constraint_num; ++ i) {
    for (int j = 0; j < _variable_num; ++ j) {
      std::cout << _coefficients(i, j) << " " ;
    }
    std::cout << " <= " << _constants[i] << std::endl ;
  }
}

/*******************************************************************************
 * Prints the {idx}th constraint to the terminal
*******************************************************************************/
void Polyhedron::PrintConstraint(int idx) const {
  std::cout << get_coefficients().row(idx) << " "
    << GetConstant(idx) << std::endl ;
}

/*******************************************************************************
 * Set a ray for the corresponding constraint. 
 * @para i the index of the constraint.
 * @para rayDirect the ray direction
*******************************************************************************/
/*
void Polyhedron::SetWitnessRay(int idx, const Vector& rayDirect) {
  _witness_ray.col(idx) = rayDirect ;
}
Vector Polyhedron::GetWitnessRay(int idx) const {
  return _witness_ray.col(idx) ; 
}
*/
/*******************************************************************************
 * Compute the witness points if it has not been computed, otherwise return
 * the witness points.
 * @para simple if simple is false, the witness point should only satisfy its
 * corresponding constraint and violate the others.
*******************************************************************************/
std::vector<Point> Polyhedron::GetWitness() {
  if ( ! _witness_point.empty() ) {
    return _witness_point ;
  } 
  if ( ! IsMinimized() ) {
    Minimize() ;
  }
  std::vector<int> activeVec = GetActiveIdx() ; 
  int consNum = activeVec.size() ;
  std::vector<Point> witness ;  
  for (int i = 0; i < consNum; ++ i) {
    witness.push_back( _witness_point[ activeVec[i] ] ) ;
  }
  return witness ;
}

/*******************************************************************************
 * Computes the exact rational solution with Flint
 * @return true if the LP has optimal solution
 * TODO These functions need to be changed (remove the output and provide functions
 * to show the result...)
*******************************************************************************/
bool Polyhedron::GetExactSolution() {
  int colNum = get_variable_num() ;
  Vector obj(colNum) ;
  // Eigen random matrix doesn't work well, so use std
  std::random_device rd ;
  std::mt19937_64 gen( rd() ) ;
  std::uniform_int_distribution<> dis(-50, 50) ;
  for (int j = 0; j < colNum; ++ j) {
    obj(j) = dis(gen) ;
  }
  return GetExactSolution(obj) ;
}

bool Polyhedron::GetExactSolution(const Vector& obj) {
  int rowNum = get_constraint_num() ;
  int colNum = get_variable_num() ;
  if (rowNum < colNum) {
    std::cerr << "Cannot solve linear equation (multiple solutions)." << std::endl ;
  }
  GlpkInterface glpkInter ;
  // minimize the obj
  bool glpkRes = glpkInter.Simplex(*this, obj) ;
  if (glpkRes == false) {
    std::cout << "Polyhedron::GetExactSolution(): LP has no optimal solution." << std::endl << std::endl ;
    return false ;
  }
  glpkInter.GetBasis() ;
  // get constraints matrix from eigen
  // The matrix is:
  //       col1  col2  ...  coln - row1    row2  ...  rown <= constant
  // row1  coef  coef       coef    -1       0           0         c   
  // row2  coef  coef       coef    0       -1           0         c
  // rown  coef  coef       coef    0       0           -1         c
  RMatrix matrix(rowNum, colNum+rowNum+1) ;
  int tmp ;
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum; ++ j) {
      tmp = Tool::GetIntOfFloat( GetCoef(i, j) ) ;
      matrix.at(i, j) = std::move(tmp) ;
    }
  }
  for (int i = 0; i < rowNum; ++ i) {
    matrix.at(i, colNum+i) = -1 ;
    tmp = Tool::GetIntOfFloat( GetConstant(i) ) ;
    matrix.at(i, colNum+rowNum) = std::move(tmp) ;
  }
  // set objective function
  RMatrix objMatrix(1, colNum+rowNum+1) ; 
  for (int j = 0; j < colNum+rowNum+1; ++ j) {
    if (j < colNum) {
      tmp = Tool::GetIntOfFloat( obj(j) ) ;
      objMatrix.at(0, j) = tmp ;
    }
    else {
      objMatrix.at(0, j) = 0 ;
    }
  }
  // solve linear equation
  std::vector<int> basic = glpkInter.get_basic_idx() ; 
  RMatrix coeff(rowNum, rowNum) ;
  RMatrix constant(rowNum, 1) ;
  for (int i = 0; i < rowNum; ++ i) {
    for(int j = 0; j < rowNum; ++ j) {
      coeff.at(i, j) = matrix.at( i, basic[j] ) ;
    }
    constant.at(i, 0) = matrix.at(i, rowNum+colNum) ;
  }
  RMatrix exactPoint( coeff.solve_dixon(constant) ) ; 
  RNumber exactRes ;
  exactRes.integer(0) ;
  for (int i = 0; i < (int)basic.size(); ++ i) {
    if (basic[i] < colNum) {
      exactRes += exactPoint.at(i, 0) * objMatrix.at(0, basic[i]) ;
    }
  }
  std::cout << "Glpk solution is: " << glpkInter.get_simplex_optimal() << std::endl ; 
  std::cout << "The exact solution is: " << exactRes << std::endl ;
  // check the result
  RMatrix objCoeff((int)basic.size(), 1) ;
  for (int i = 0; i < (int)basic.size(); ++ i) {
    objCoeff.at(i, 0) = objMatrix.at(0, basic[i]) ;
  }
  RMatrix lamda( coeff.transpose().solve_dixon(objCoeff) );  
  //RMatrix newObj(lamda.transpose() * matrix - objMatrix) ; 
  RMatrix newObj(objMatrix - lamda.transpose() * matrix) ; 
  std::vector<int> nonbasic = glpkInter.get_non_basic_idx() ; 
  RNumber zero ;
  zero.set_zero() ;
  int currIdx, state ;
  for (int i = 0; i < (int) nonbasic.size(); ++ i) { 
    currIdx = nonbasic[i] ;
    state = glpkInter.GetVariState(currIdx) ;
    // we only use minimize in our program
    // max obj is min - obj
    /*
    if (objdir == GLP_MAX) {
      if (state == GLP_NF && newObj.at(0, currIdx) == zero) {
        continue ;
      }
      else if (state == GLP_NL && newObj.at(0, currIdx) <= zero) {
        continue ;
      }
      else if (state == GLP_NU && newObj.at(0, currIdx) >= zero) {
        continue ;
      }
      else {
        std::cerr << "Incorrect simplex result." << std::endl ; 
        return false ;
      }
    }
    else {*/
      if (state == GLP_NF && newObj.at(0, currIdx) == zero) {
        continue ;
      }
      else if (state == GLP_NL && newObj.at(0, currIdx) >= zero) {
        continue ;
      }
      else if (state == GLP_NU && newObj.at(0, currIdx) <= zero) {
        continue ;
      }
      else {
        std::cerr << "Incorrect simplex result." << std::endl ; 
        return false ;
      }
    //}
  }
  std::cout << "The simplex result is correct." << std::endl ;
  return true ; 
}
