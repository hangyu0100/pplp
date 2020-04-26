/*******************************************************************************
 * Copyright (c) 2017 Jan. Verimag. All rights reserved.
 * @author Hang YU
*******************************************************************************/

#include <string>
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include "pplp_apronInterface.h"

using namespace PPLP ;

/*******************************************************************************
 * See apronInterface.h
*******************************************************************************/
ApronInterface::ApronInterface(const Polyhedron& poly) {
  _man = pk_manager_alloc(true) ;
  _vari_num = poly.get_variable_num() ;
  int consNum = poly.get_constraint_num() ;
  _name_of_vari = new char*[_vari_num] ; 
  std::string variName ;
  for (int j = 0; j < _vari_num; ++ j) {
#ifndef NO_STRING
    variName = "x" + std::to_string(j) ;
#endif
    _name_of_vari[j] = new char[ variName.length()+1 ] ;
    std::strcpy( _name_of_vari[j], variName.c_str() ) ;
  }
  _env = ap_environment_alloc(NULL, 0, (ap_var_t*)_name_of_vari, _vari_num) ;
  int val, constant ;
  ap_linexpr1_t expr ;
  ap_lincons1_t cons ;
  // set constraints
  _cons_array = ap_lincons1_array_make(_env, consNum) ; 
  for (int i = 0; i < consNum; ++ i) {
    expr = ap_linexpr1_make(_env, AP_LINEXPR_DENSE, 0) ;
    for (int j = 0; j < _vari_num; ++ j) {
      val = Tool::GetIntOfFloat( -poly.GetCoef(i, j) ) ;
      ap_linexpr1_set_coeff_scalar_double(&expr, _name_of_vari[j], val) ;
    }
    constant = Tool::GetIntOfFloat( poly.GetConstant(i) ) ;
    ap_linexpr1_set_cst_scalar_double(&expr, constant) ;
    cons = ap_lincons1_make(AP_CONS_SUPEQ, &expr, NULL) ;
    ap_lincons1_array_set(&_cons_array, i, &cons) ;
  }
  _abs = ap_abstract1_of_lincons_array(_man, _env, &_cons_array) ;
  _generator = ap_abstract1_to_generator_array(_man, &_abs) ;
}

ApronInterface::ApronInterface(const RMatrix& optimals) {
  _man = pk_manager_alloc(true) ;
  // each row is a constraint, in form of AX + b >= 0
  int consNum = optimals.rows() ;
  // the last column is the constant
  _vari_num = optimals.cols() - 1 ;
  _name_of_vari = new char*[_vari_num] ; 
  std::string variName ;
  for (int j = 0; j < _vari_num; ++ j) {
#ifndef NO_STRING
    variName = "x" + std::to_string(j) ;
#endif
    _name_of_vari[j] = new char[ variName.length()+1 ] ;
    std::strcpy( _name_of_vari[j], variName.c_str() ) ;
  }
  _env = ap_environment_alloc(NULL, 0, (ap_var_t*)_name_of_vari, _vari_num) ;
  RNumber val ;
  ap_linexpr1_t expr ;
  ap_lincons1_t cons ;
  // set constraints
  _cons_array = ap_lincons1_array_make(_env, consNum) ; 
  ap_scalar_t* currScalar = ap_scalar_alloc() ;
  for (int i = 0; i < consNum; ++ i) {
    expr = ap_linexpr1_make(_env, AP_LINEXPR_DENSE, 0) ;
    for (int j = 0; j < _vari_num; ++ j) {
      val = optimals.at(i, j) ;
      FlintToApron(currScalar, val) ;
      ap_linexpr1_set_coeff_scalar(&expr, _name_of_vari[j], currScalar) ;
    }
    val = optimals.at(i, _vari_num) ;
    FlintToApron(currScalar, val) ;
    ap_linexpr1_set_cst_scalar(&expr, currScalar) ;
    cons = ap_lincons1_make(AP_CONS_SUPEQ, &expr, NULL) ;
    ap_lincons1_array_set(&_cons_array, i, &cons) ;
  }
  ap_scalar_free(currScalar) ;
  _abs = ap_abstract1_of_lincons_array(_man, _env, &_cons_array) ;
  _generator = ap_abstract1_to_generator_array(_man, &_abs) ;
}

ApronInterface::~ApronInterface() {
  ap_lincons1_array_clear(&_cons_array) ;
  ap_generator1_array_clear(&_generator) ;
  ap_abstract1_clear(_man, &_abs) ;
  ap_environment_free(_env) ;
  ap_manager_free(_man) ;
  for (int i = 0; i < _vari_num; ++ i) {
    delete[] _name_of_vari[i] ;
  }
  delete[] _name_of_vari ;
}

/*******************************************************************************
 * @return the polyhedron which is obtained from the apron abstract domain
*******************************************************************************/
Polyhedron ApronInterface::GetPoly () {
  return GetPoly(_abs) ;
}

Polyhedron ApronInterface::GetPoly (ap_abstract1_t abs) {
  ap_lincons1_array_t consArray = ap_abstract1_to_lincons_array(_man, &abs) ;
  int consNum = ap_lincons1_array_size(&consArray) ; 
  ap_scalar_t* currScalar = ap_scalar_alloc() ;
  ap_coeff_t* currVal = ap_coeff_alloc_set_scalar(currScalar) ;
  ap_linexpr1_t expr = ap_linexpr1_make(_env, AP_LINEXPR_DENSE, 0) ;
  ap_lincons1_t cons = ap_lincons1_make(AP_CONS_SUPEQ, &expr, NULL) ;
  double currCoeff ;
  Polyhedron newPoly(consNum, _vari_num) ;
  ap_constyp_t op ;
  for (int i = 0; i < consNum; ++ i) {
    cons = ap_lincons1_array_get(&consArray, i) ;
    for (int j = 0; j < _vari_num; ++ j) {
      ap_lincons1_get_coeff (currVal, &cons, _name_of_vari[j]) ;
      ap_double_set_scalar(&currCoeff, currVal->val.scalar, GMP_RNDN) ;
      newPoly.SetCoef(i, j, - currCoeff) ;
    }
    ap_lincons1_get_cst(currVal, &cons) ; 
    ap_double_set_scalar(&currCoeff, currVal->val.scalar, GMP_RNDN) ;
    newPoly.SetConstant(i, currCoeff) ;
    op = *ap_lincons1_constypref(&cons) ;
    if (op == AP_CONS_SUPEQ) {
      //newPoly.SetOperator(i, ConstraintOperator::lesseq) ;
      //todo
    }
    else if (op == AP_CONS_EQ) {
      //newPoly.SetOperator(i, ConstraintOperator::equal) ;
    }
    /*
    // if there are other operators
    else if () ...
     */
  }
  return newPoly ;
}

/*******************************************************************************
 * @return true if the two polyhedra contain the same constraints
*******************************************************************************/
bool ApronInterface::CompareAbstract(const ApronInterface& apron) {
  ap_abstract1_t apronAbs = apron.get_abs() ;
  return ap_abstract1_is_eq( _man, &_abs, &apronAbs) ;  
}

/*******************************************************************************
 * @return the number of generators
*******************************************************************************/
int ApronInterface::GetGeneratorNum() {
  int size = (int)ap_generator1_array_size(&_generator) ;
  return size ;
}

ap_abstract1_t ApronInterface::get_abs() const {
  return _abs ;
}


/*******************************************************************************
 * Print the irredundant constraints
*******************************************************************************/
void ApronInterface::Print() {
  ap_abstract1_fprint(stdout, _man, &_abs) ;
}

/*******************************************************************************
 * @return the generators as a 1-D array
*******************************************************************************/
std::vector<double> ApronInterface::get_gen_array() {
  ap_generator1_t currGen ;
  ap_linexpr1_t expr ;
  ap_scalar_t* currScalar = ap_scalar_alloc() ;
  ap_coeff_t* currVal = ap_coeff_alloc_set_scalar(currScalar) ;
  double currCoeff ;
  int size = ap_generator1_array_size(&_generator) ;
  for (int i = 0; i < size; ++ i) {
    currGen = ap_generator1_array_get(&_generator, i) ;
    expr = ap_generator1_linexpr1ref(&currGen) ;
    for (int j = 0; j < _vari_num; ++ j) {
      currVal = ap_linexpr1_coeffref (&expr, _name_of_vari[j]) ;
      ap_double_set_scalar(&currCoeff, currVal->val.scalar, GMP_RNDN) ;
      _gen_array.push_back(currCoeff) ;
    }
  }
  return _gen_array ;
}

void ApronInterface::PrintGenerator() {
  ap_generator1_array_fprint(stdout, &_generator) ;
}

/*******************************************************************************
 * @para idx the index of the generator
 * @para dim the index of the dimension
 * @return the value of {idx}th generator in {dim} dimension 
*******************************************************************************/
double ApronInterface::GetGenVal(int idx, int dim) {
  return _gen_array[idx * _vari_num + dim] ;
}

/*******************************************************************************
 * Projects the first num variables
 * @para idx the index of variables to be projected
 * @para matrix the projected polyhedron matrix to be compared
*******************************************************************************/
bool ApronInterface::CmpProjectedPoly(const std::vector<int>& idx, 
    const RMatrix& matrix) {
  int num = idx.size() ;
  char** names = new char*[num] ; 
  std::string variName ;
  for (int j = 0; j < num; ++ j) {
#ifndef NO_STRING
    variName = "x" + std::to_string( idx[j] ) ;
#endif
    names[j] = new char[ variName.length()+1 ] ;
    std::strcpy( names[j], variName.c_str() ) ;
  }

  ap_abstract1_t projAbs = 
      ap_abstract1_forget_array(_man, false, &_abs, (ap_var_t*)names, num, false) ;

  for (int i = 0; i < num; ++ i) {
    delete[] names[i] ;
  }
  delete[] names ;

  ApronInterface currApron(matrix) ;
  //ap_abstract1_fprint(stdout, _man, &projAbs) ;
  //ap_abstract1_fprint(stdout, _man, &currApron._abs) ;
  bool res = ap_abstract1_is_eq(_man, &projAbs, &currApron._abs) ;
  // for test
  if (res == false) {
    ap_abstract1_fprint(stdout, _man, &projAbs) ;
  }
  return res ;
}

void ApronInterface::Project(int projNum) { 
  char** names = new char*[projNum] ; 
  std::string variName ;
  for (int j = 0; j < projNum; ++ j) {
    variName = "x" + std::to_string(j) ;
    names[j] = new char[ variName.length()+1 ] ;
    std::strcpy( names[j], variName.c_str() ) ;
  }

  ap_abstract1_t projAbs = 
      ap_abstract1_forget_array(_man, false, &_abs, (ap_var_t*)names, projNum, false) ;

  #ifdef APRON_ASK_CONS
    ap_abstract1_to_lincons_array(_man, &projAbs) ;
  #endif

  for (int i = 0; i < projNum; ++ i) {
    delete[] names[i] ;
  }
  delete[] names ;

  ap_abstract1_fprint(stdout, _man, &projAbs) ;
}

bool ApronInterface::CmpConvexHull(const Polyhedron& poly2, const RMatrix& matrix) {
  int val, constant ;
  ap_linexpr1_t expr ;
  ap_lincons1_t cons ;
  int consNum2 = poly2.get_constraint_num() ;
  ap_lincons1_array_t _cons_array2 = ap_lincons1_array_make(_env, consNum2) ; 
  for (int i = 0; i < consNum2; ++ i) {
    expr = ap_linexpr1_make(_env, AP_LINEXPR_DENSE, 0) ;
    for (int j = 0; j < _vari_num; ++ j) {
      val = Tool::GetIntOfFloat( -poly2.GetCoef(i, j) ) ;
      ap_linexpr1_set_coeff_scalar_double(&expr, _name_of_vari[j], val) ;
    }
    constant = Tool::GetIntOfFloat( poly2.GetConstant(i) ) ;
    ap_linexpr1_set_cst_scalar_double(&expr, constant) ;
    cons = ap_lincons1_make(AP_CONS_SUPEQ, &expr, NULL) ;
    ap_lincons1_array_set(&_cons_array2, i, &cons) ;
  }
  ap_abstract1_t abs2 = ap_abstract1_of_lincons_array(_man, _env, &_cons_array2) ;
  
  ap_abstract1_t resAbs = ap_abstract1_join(_man , false, &_abs, &abs2) ;
  
  ApronInterface currApron(matrix) ;
  bool res = ap_abstract1_is_eq(_man, &resAbs, &currApron._abs) ;
  if (res == false) {
    ap_abstract1_fprint(stdout, _man, &resAbs) ;
  }
  return res ;
}

ap_lincons1_array_t ApronInterface::ConvexHull(const Polyhedron& poly2) {
  int val, constant ;
  ap_linexpr1_t expr ;
  ap_lincons1_t cons ;
  int consNum2 = poly2.get_constraint_num() ;
  ap_lincons1_array_t _cons_array2 = ap_lincons1_array_make(_env, consNum2) ; 
  for (int i = 0; i < consNum2; ++ i) {
    expr = ap_linexpr1_make(_env, AP_LINEXPR_DENSE, 0) ;
    for (int j = 0; j < _vari_num; ++ j) {
      val = Tool::GetIntOfFloat( -poly2.GetCoef(i, j) ) ;
      ap_linexpr1_set_coeff_scalar_double(&expr, _name_of_vari[j], val) ;
    }
    constant = Tool::GetIntOfFloat( poly2.GetConstant(i) ) ;
    ap_linexpr1_set_cst_scalar_double(&expr, constant) ;
    cons = ap_lincons1_make(AP_CONS_SUPEQ, &expr, NULL) ;
    ap_lincons1_array_set(&_cons_array2, i, &cons) ;
  }
  ap_abstract1_t abs2 = ap_abstract1_of_lincons_array(_man, _env, &_cons_array2) ;
  ap_abstract1_t resAbs = ap_abstract1_join(_man , false, &_abs, &abs2) ;
  ap_lincons1_array_t res = ap_abstract1_to_lincons_array(_man, &resAbs) ;

  ap_generator1_array_t generator = ap_abstract1_to_generator_array(_man, &resAbs) ;
  int genSize = (int)ap_generator1_array_size(&generator) ;
  std::cout << "generator," << genSize << "," ;

  return res ;

}

/*******************************************************************************
 * Convert RNumber to ap_scalar_t
*******************************************************************************/
void ApronInterface::FlintToApron(ap_scalar_t* scalar, const RNumber& val) {
  mpq_t currR ;
  mpq_init(currR) ;
  mpq_set_str( currR, val.to_string().c_str(), 10 ) ;
  ap_scalar_set_mpq(scalar, currR) ;
  mpq_clear(currR) ;
}

int ApronInterface::GetConsNum() {
  ap_lincons1_array_t arrayCons = ap_abstract1_to_lincons_array(_man, &_abs) ;
  return ap_lincons1_array_size(&arrayCons) ;
}
