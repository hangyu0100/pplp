#include <iostream>
#include <random>
#include "pplp_polyhedron.h"
#include "pplp_ioInterface.h"
#include "pplp_tbbParallel.h"

using namespace PPLP ;

bool TestConstructor() {
  int consNum = 3 ;
  int variNum = 3 ;
  int eqConsNum = 1 ;
  Matrix ineqCoeff(consNum, variNum+1) ;
  ineqCoeff << 3, 0, -1, 12,
           -1, 0, -3, 7,
           -1, 0, 2, -4 ;
  Matrix eqCoeff(eqConsNum, variNum+1) ;
  eqCoeff << 0, 1, 0, 2 ;

  RMatrix ineqCoeffR(consNum, variNum+1) ;
  RMatrix eqCoeffR(eqConsNum, variNum+1) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum+1; ++ j) {
      ineqCoeffR.at(i,j) = (int)round(ineqCoeff(i,j)) ;
    }
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    for (int j = 0; j < variNum+1; ++ j) {
      eqCoeffR.at(i,j) = (int)eqCoeff(i,j) ;
    }
  } 

  Polyhedron poly1 ;
  assert(poly1.get_constraint_num() == 0) ;
  assert(poly1.get_eq_constraint_num() == 0) ;
  assert(poly1.get_variable_num() == 0) ;
  assert( poly1.IsTop() ) ;
  poly1.SetSize(consNum, variNum, eqConsNum) ;
  assert(poly1.get_constraint_num() == 3) ;
  assert(poly1.get_eq_constraint_num() == 1) ;
  assert(poly1.get_variable_num() == 3) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      poly1.SetCoef( i, j, ineqCoeff(i,j) ) ;
    }
    poly1.SetConstant( i, ineqCoeff(i,variNum) ) ;
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      poly1.SetEqCoef( i, j, eqCoeff(i,j) ) ;
    }
    poly1.SetEqConstant( i, eqCoeff(i,variNum) ) ;
  }

  Polyhedron poly4(consNum, variNum, eqConsNum) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      poly4.SetCoef( i, j, ineqCoeff(i,j) ) ;
    }
    poly4.SetConstant( i, ineqCoeff(i,variNum) ) ;
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      poly4.SetEqCoef( i, j, eqCoeff(i,j) ) ;
    }
    poly4.SetEqConstant( i, eqCoeff(i,variNum) ) ;
  }
  assert( poly1.get_coefficients() == poly4.get_coefficients() && poly1.get_constants() == poly4.get_constants() ) ;
  assert( poly1.get_eq_coefficients() == poly4.get_eq_coefficients() && poly1.get_eq_constants() == poly4.get_eq_constants() ) ;
  assert( poly1.get_coefficients() == poly4.get_coefficients() && poly1.get_constants() == poly4.get_constants() ) ;
  assert( poly1.get_eq_coefficients() == poly4.get_eq_coefficients() && poly1.get_eq_constants() == poly4.get_eq_constants() ) ;

  Polyhedron poly2(3, 3, 1) ; 
  assert(poly2.get_constraint_num() == 3) ;
  assert(poly2.get_eq_constraint_num() == 1) ;
  assert(poly2.get_variable_num() == 3) ;
  poly2.set_coefficients( ineqCoeff.block(0, 0, consNum, variNum) ) ;
  poly2.set_constants( ineqCoeff.col(variNum).transpose() ) ;
  poly2.set_eq_coefficients( eqCoeff.block(0, 0, eqConsNum, variNum) ) ;
  poly2.set_eq_constants( eqCoeff.col(variNum).transpose() ) ;

  Polyhedron poly3( ineqCoeffR, std::vector<int>(), eqCoeffR, std::vector<int>() ) ;

  assert( poly1.get_coefficients() == poly2.get_coefficients() && poly1.get_constants() == poly2.get_constants() ) ;
  assert( poly1.get_eq_coefficients() == poly2.get_eq_coefficients() && poly1.get_eq_constants() == poly2.get_eq_constants() ) ;
  assert( poly1.get_coefficients() == poly3.get_coefficients() && poly1.get_constants() == poly3.get_constants() ) ;
  assert( poly1.get_eq_coefficients() == poly3.get_eq_coefficients() && poly1.get_eq_constants() == poly3.get_eq_constants() ) ;

  poly2.Init() ;
  poly2.set_coefficients( ineqCoeff.block(0, 0, consNum, variNum) ) ;
  poly2.set_constants( ineqCoeff.col(variNum).transpose() ) ;
  poly2.set_eq_coefficients( eqCoeff.block(0, 0, eqConsNum, variNum) ) ;
  poly2.set_eq_constants( eqCoeff.col(variNum).transpose() ) ;
  assert( poly1.get_coefficients() == poly2.get_coefficients() && poly1.get_constants() == poly2.get_constants() ) ;
  assert( poly1.get_eq_coefficients() == poly2.get_eq_coefficients() && poly1.get_eq_constants() == poly2.get_eq_constants() ) ;

  return true ;
}

bool TestEquality() {
  Polyhedron poly1(2, 3, 1) ;
  Polyhedron poly2(2, 3, 1) ;

  poly1.SetCoef(0, 1, 2) ;
  poly1.SetCoef(0, 2, 1) ;
  poly1.SetConstant(0, -1) ;
  poly2.SetCoef(1, 1, 2) ;
  poly2.SetCoef(1, 2, 1) ;
  poly2.SetConstant(1, -1) ;
  poly1.SetCoef(1, 0, -1) ;
  poly1.SetCoef(1, 2, 3) ;
  poly1.SetConstant(1, 2) ;
  poly2.SetCoef(0, 0, -1) ;
  poly2.SetCoef(0, 2, 3) ;
  poly2.SetConstant(0, 2) ;

  poly1.SetEqCoef(0, 0, 1) ;
  poly1.SetEqCoef(0, 1, 2) ;
  poly1.SetEqCoef(0, 2, -2) ;
  poly1.SetEqConstant(0, -1) ;
  poly2.SetEqCoef(0, 0, 1) ;
  poly2.SetEqCoef(0, 1, 2) ;
  poly2.SetEqCoef(0, 2, -2) ;
  poly2.SetEqConstant(0, -1) ;

  return poly1.IsEqualTo(poly2) ;
}

bool TestSubPoly() {
  int consNum = 4 ;
  int variNum = 3 ;
  int eqConsNum = 4 ;
  Polyhedron poly(consNum, variNum, eqConsNum) ;
  srand( time(NULL) ) ;

  std::random_device rd ;
  std::mt19937 gen(rd()) ;
  std::uniform_int_distribution<> dis(-10, 10) ;

  int coeff ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      coeff = dis(gen) ;
      poly.SetCoef(i, j, coeff) ;
    }
    coeff = dis(gen) ;
    poly.SetConstant(i, coeff) ;
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      coeff = dis(gen) ;
      poly.SetEqCoef(i, j, coeff) ;
    }
    coeff = dis(gen) ;
    poly.SetEqConstant(i, coeff) ;
  }

  std::vector<int> idx1 ;
  idx1.push_back(1) ;
  idx1.push_back(2) ;
  Polyhedron sub1 = poly.GetSubPoly(idx1) ;

  assert( sub1.get_constraint_num() == idx1.size() ) ;
  for (unsigned i = 0; i < idx1.size(); ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub1.GetCoef(i, j) == poly.GetCoef(idx1[i], j) ) ;
    }
    assert( sub1.GetConstant(i) == poly.GetConstant( idx1[i] ) ) ;
  }
  for (int i = 0; i < eqConsNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub1.GetEqCoef(i, j) == poly.GetEqCoef(i, j) ) ;
    }
    assert( sub1.GetEqConstant(i) == poly.GetEqConstant(i) ) ;
  }

  std::vector<int> idx2 ;
  idx2.push_back(0) ;
  idx2.push_back(3) ;
  Polyhedron sub2 = poly.GetSubPoly(std::vector<int>(), idx2) ;
  assert( sub2.get_eq_constraint_num() == idx2.size() ) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub2.GetCoef(i, j) == poly.GetCoef(i, j) ) ;
    }
    assert( sub2.GetConstant(i) == poly.GetConstant(i) ) ;
  }
  for (unsigned i = 0; i < idx2.size(); ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub2.GetEqCoef(i, j) == poly.GetEqCoef(idx2[i], j) ) ;
    }
    assert( sub2.GetEqConstant(i) == poly.GetEqConstant( idx2[i] ) ) ;
  }

  Polyhedron sub3 = poly.GetSubPoly(idx1, idx2) ;
  assert( sub3.get_constraint_num() == idx1.size() ) ;
  assert( sub3.get_eq_constraint_num() == idx2.size() ) ;
  for (unsigned i = 0; i < idx1.size(); ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub3.GetCoef(i, j) == poly.GetCoef(idx1[i], j) ) ;
    }
    assert( sub3.GetConstant(i) == poly.GetConstant( idx1[i] ) ) ;
  }
  for (unsigned i = 0; i < idx2.size(); ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      assert( sub3.GetEqCoef(i, j) == poly.GetEqCoef(idx2[i], j) ) ;
    }
    assert( sub3.GetEqConstant(i) == poly.GetEqConstant( idx2[i] ) ) ;
  }
  
  return true ;
}

bool TestIo() {
  int consNum = 3 ;
  int variNum = 3 ;
  int eqConsNum = 1 ;
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra("./testfiles/test.poly") ;
  assert( ioInter.get_cons_num() == consNum+eqConsNum) ;
  assert( ioInter.get_vari_num() == variNum) ;
  Polyhedron poly = polyVec[0] ;
  Matrix ineqCoeff(consNum, variNum+1) ;
  ineqCoeff << 3, 0, -1, 12,
           -1, 0, -3, 7,
           -1, 0, 2, -4 ;
  Matrix eqCoeff(eqConsNum, variNum+1) ;
  eqCoeff << 0, 1, 0, 2 ;

  assert( poly.get_coefficients() == ineqCoeff.block(0, 0, consNum, variNum) ) ;
  assert( poly.get_constants() == ineqCoeff.col(variNum).transpose() ) ;
  assert( poly.get_eq_coefficients() == eqCoeff.block(0, 0, eqConsNum, variNum) ) ;
  assert( poly.get_eq_constants() == eqCoeff.col(variNum).transpose() ) ;

  return true ;
}

bool TestGaussian() {
  RMatrix matrix(4, 5) ;
  matrix.at(0,0) = 2 ;
  matrix.at(0,1) = 3 ;
  matrix.at(0,2) = -2 ;
  matrix.at(0,3) = 1 ;
  matrix.at(0,4) = 2 ;
  matrix.at(1,0) = 3 ;
  matrix.at(1,1) = -2 ;
  matrix.at(1,2) = 4 ;
  matrix.at(1,3) = -1 ;
  matrix.at(1,4) = 3 ;
  matrix.at(2,0) = 10 ;
  matrix.at(2,1) = 15 ;
  matrix.at(2,2) = -10 ;
  matrix.at(2,3) = 5 ;
  matrix.at(2,4) = 10 ;
  matrix.at(3,0) = -1 ;
  matrix.at(3,1) = -4 ;
  matrix.at(3,2) = 2 ;
  matrix.at(3,3) = 3 ;
  matrix.at(3,4) = 1 ;

  ReducedMatrix reduced = Tool::GaussianElimination(matrix) ;
  assert(reduced.succeed) ;
  int consNum = reduced.rowIdx.size() ;
  int variNum = reduced.colIdx.size() ;
  for (int j = 0; j < variNum+1; ++ j) {
    assert( reduced.matrix.at(2,j).is_zero() ) ;
  }
  int rowIdx ;
  RMatrix rm(consNum, variNum+1) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum+1; ++ j) {
      rowIdx = reduced.rowIdx[i] ;
      rm.at(i,j) = reduced.matrix.at(rowIdx, j) ;
    }
  }
  for (int i = 0; i < consNum; ++ i) {
    assert( rm.at(i,i).is_one() ) ;
  }

  matrix.at(2,4) = 9 ;
  ReducedMatrix reduced2 = Tool::GaussianElimination(matrix) ;
  assert( ! reduced2.succeed) ;

  return true ;
}


bool TestGaussianFloat() {
  Matrix matrix(4, 5) ;
  matrix << 2, 3, -2, 1, 2,
            3, -2, 4, -1, 3,
            10, 15, -10, 5, 10,
            -1, -4, 2, 3, 1 ;

  ReducedMatrixFloat reduced = Tool::GaussianEliminationFloat(matrix) ;
  assert(reduced.succeed) ;
  int consNum = reduced.rowIdx.size() ;
  int variNum = reduced.colIdx.size() ;
  for (int j = 0; j < variNum+1; ++ j) {
    assert(reduced.matrix(2,j) == 0) ;
  }
  int rowIdx ;
  Matrix rm(consNum, variNum+1) ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum+1; ++ j) {
      rowIdx = reduced.rowIdx[i] ;
      rm(i,j) = reduced.matrix(rowIdx, j) ;
    }
  }
  for (int i = 0; i < consNum; ++ i) {
    assert(rm(i,i) == 1) ;
  }


  matrix(2,4) = 9 ;
  ReducedMatrixFloat reduced2 = Tool::GaussianEliminationFloat(matrix) ;
  assert( ! reduced2.succeed) ;

  matrix(2,4) = 10 ;
  matrix.row(2) = matrix.row(2) / 3 ;
  ReducedMatrixFloat reduced3 = Tool::GaussianEliminationFloat(matrix) ;
  assert(reduced3.succeed) ;

  return true ;
}

bool TestSubsEq() {
  int consNum = 3 ;
  int variNum = 3 ;
  int eqConsNum = 3 ;
  Matrix ineqCoeff(consNum, variNum+1) ;
  ineqCoeff << 3, 2, -1, 12,
           -1, 1, -3, 7,
           -3, -1, 1, -7 ;
  Matrix eqCoeff(eqConsNum, variNum+1) ;
  eqCoeff << 0, 1, 0, 2,
            0, 2, 3, -2,
            0, 6, 9, -6 ;
  Polyhedron poly(consNum, variNum, eqConsNum) ;
  poly.set_coefficients( ineqCoeff.block(0, 0, consNum, variNum) ) ;
  poly.set_constants( ineqCoeff.col(variNum).transpose() ) ;
  poly.set_eq_coefficients( eqCoeff.block(0, 0, eqConsNum, variNum) ) ;
  poly.set_eq_constants( eqCoeff.col(variNum).transpose() ) ;
  Polyhedron poly2 = poly ;

  ReducedMatrix redMatrixR = poly.SubsEqConsR() ;
  ReducedMatrixFloat redMatrixF = poly2.SubsEqConsF() ;
  assert(poly.GetEqCoef(0,1) == 1) ;
  assert(poly.GetEqCoef(1,2) == 1) ;
  assert( poly.get_eq_coefficients() == poly2.get_eq_coefficients() ) ;
  assert( poly.get_eq_constants() == poly2.get_eq_constants() ) ;
  assert( poly.get_eq_constraint_num() == poly2.get_eq_constraint_num() ) ;

  poly.SubsIneqConsR(redMatrixR) ;
  poly2.SubsIneqConsF(redMatrixF) ;
  assert( poly.get_eq_coefficients() == poly2.get_eq_coefficients() ) ;
  assert( poly.get_eq_constants() == poly2.get_eq_constants() ) ;
  assert( poly.get_coefficients() == poly2.get_coefficients() ) ;
  assert( poly.get_constants() == poly2.get_constants() ) ;

  return true ;
}

bool TestMiniIneqOnly() {
  int consNum = 5 ;
  int variNum = 3 ;
  Matrix ineqCoeff(consNum, variNum+1) ;
  ineqCoeff << -21, -33, -39, 30,
              -168, -174, -67, 261,
              -273, -357, -311, 429,
              -21, -15, 10, 31,
              -168, -156, -18, 289 ;
  Polyhedron poly(consNum, variNum) ;

  poly.set_coefficients( ineqCoeff.block(0, 0, consNum, variNum) ) ;
  poly.set_constants( ineqCoeff.col(variNum).transpose() ) ;

  poly.Minimize() ;
  std::vector<int> activeIdx = poly.GetActiveIdx() ;
  assert(activeIdx[0] == 0) ;
  assert(activeIdx[1] == 3) ;
  
  return true ;
}

bool TestMinimize() {
  int consNum = 5 ;
  int eqConsNum = 2 ;
  int variNum = 5 ;
  Matrix ineqCoeff(consNum, variNum+1) ;
  ineqCoeff << 12, -11, -21, -33, -39, 30,
              21, -3, -168, -174, -67, 261,
              -5, 10, -273, -357, -311, 429,
              -21, 10, -21, -15, 10, 31,
              -3, 21, -168, -156, -18, 289 ;
  Polyhedron poly(consNum, variNum, eqConsNum) ;
  Matrix eqCoeff(eqConsNum, variNum+1) ;
  eqCoeff << 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0 ;

  poly.set_coefficients( ineqCoeff.block(0, 0, consNum, variNum) ) ;
  poly.set_constants( ineqCoeff.col(variNum).transpose() ) ;
  poly.set_eq_coefficients( eqCoeff.block(0, 0, eqConsNum, variNum) ) ;
  poly.set_eq_constants( eqCoeff.col(variNum).transpose() ) ;

  poly.Minimize() ;
  std::vector<int> activeIdx = poly.GetActiveIdx() ;
  assert(activeIdx[0] == 0) ;
  assert(activeIdx[1] == 3) ;
  
  return true ;
}

bool TestConstructorProj() {
  int consNum = 3, variNum = 2 ;
  Polyhedron poly(consNum, variNum) ;
  Matrix coeff(consNum, variNum) ;
  Vector constant(consNum) ;
  coeff << -3, 1,
          -1, -2, 
          4, 1 ;
  constant << -1, -5, 13 ;
  poly.set_coefficients( std::move(coeff) ) ;
  poly.set_constants( std::move(constant) ) ;

  TbbParallel tbb(poly) ;
  // the second parameter 0 enforces the algo choose a start task point
  // this is for test
  tbb.PlpParallel(1, 0) ;
  // the result should be x>=1, x<=5
  RMatrix res( tbb.GetOptimalMatrix() ) ; 
  RNumber ratio1( -res.at(0,2) / res.at(0,1) ) ; 
  RNumber ratio2( -res.at(1,2) / res.at(1,1) ) ; 
  RNumber five, zero ;
  five.set_integer(5) ;
  zero.set_zero() ;

  if ( ratio1.is_one() ) {
    assert(ratio2 == five) ;
    assert(res.at(0,1) > zero) ;
    assert(res.at(0,2) < zero) ;
    assert(res.at(1,1) < zero) ;
    assert(res.at(1,2) > zero) ;
  }
  else if ( ratio2.is_one() ) {
    assert(ratio1 == five) ;
    assert(res.at(1,1) > zero) ;
    assert(res.at(1,2) < zero) ;
    assert(res.at(0,1) < zero) ;
    assert(res.at(0,2) > zero) ;
  }
  else {
    return false ;
  }

  return true ;
}

int main (int argc, char* argv[]) {

  // constructor of polyhedron
  assert( TestConstructor() ) ;
  std::cout << "Test constructor: OK" << std::endl ;

  // equality of polyhedra
  assert( TestEquality() ) ;
  std::cout << "Test polyhedra equality: OK" << std::endl ;
  
  // subpoly
  assert( TestSubPoly() ) ;
  std::cout << "Test subpoly: OK" << std::endl ;

  // io
  assert( TestIo() ) ;
  std::cout << "Test io: OK" << std::endl ;

  // minimization
  // Gaussian
  assert( TestGaussian() ) ;
  std::cout << "Test Gaussian: OK" << std::endl ;
  // Gaussian float
  assert( TestGaussianFloat() ) ;
  std::cout << "Test GaussianFloat: OK" << std::endl ;
  // Equalities
  assert( TestSubsEq() ) ;
  std::cout << "Test substitute equalities: OK" << std::endl ;
  // raytracing
  // TODO maybe refine this part
  assert( TestMiniIneqOnly() ) ;
  std::cout << "Test minimization with inequalities only: OK" << std::endl ;
  assert( TestMinimize() ) ;
  std::cout << "Test minimization with inequalities and equalities: OK" << std::endl ;

  // projection inequalities
  assert( TestConstructorProj() ) ;
  std::cout << "Test constructor of projection: OK" << std::endl ;
  
  
  
  
  return 0 ;
}
