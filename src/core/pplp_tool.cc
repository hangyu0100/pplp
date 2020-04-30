#include <limits>
#include <cmath>
#include "pplp_tool.h"

using namespace PPLP ;

#ifdef DEBUGINFO_PLP
#include <mutex>
#include <iostream>

extern std::mutex log_mtx;
#endif

double Tool::GetDotProductThreshold(const Vector& vec1, const Vector& vec2) {
  int size = vec1.size() ;
  double epsilon = std::numeric_limits<double>::epsilon() ;
  double res = 1.01 * size * epsilon * vec1.cwiseAbs() * vec2.cwiseAbs().transpose() ;
  return res ;
}

/*
double Tool::GetDotProductThreshold(const Vector& vec1, const Vector& vec2) {
  int size = vec1.size() ;
  double norm1 = vec1.norm() ;
  double norm2 = vec2.norm() ;
  return GetDotProductThreshold(size, norm1, norm2) ;
}

double Tool::GetDotProductThreshold(int size, double norm1, double norm2) {
  double epsilon = std::numeric_limits<double>::epsilon() ;
  int s = ceil( log2(size) ) + 1 ; 
  return (s * epsilon) / (1 - s * epsilon) * norm1  * norm2 ;
}
*/

/*******************************************************************************
 * Transfers rational matrix to double matrix 
 * @para rmatrix the rational matrix to be transfered
 * @return the double matrix
*******************************************************************************/
Matrix Tool::GetDoubleMatrix(const RMatrix& rmatrix) {
  int rows = rmatrix.rows() ;
  int cols = rmatrix.cols() ;
  return GetDoubleBlock(rmatrix, 0, 0, rows, cols) ;
}

Matrix Tool::GetDoubleBlock(const RMatrix& rmatrix, int beginRow, int beginCol,
      int rowNum, int colNum ) {
  if ( beginRow + rowNum > rmatrix.rows() || beginCol + colNum > rmatrix.cols() ) {
    std::cerr << "Tool::GetDoubleBlock error. Index is out of range." << std::endl ;
    std::terminate() ;
  }

  Matrix currMatrix(rowNum, colNum) ;
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum; ++ j) {
      currMatrix(i, j) = Tool::GetDouble( (RNumber)rmatrix.at(beginRow+i,beginCol+j) ) ;
    }
  } 
  return currMatrix ;
}

Vector Tool::GetDoubleVector(const RMatrix& rmatrix) {
  assert(rmatrix.cols() == 1) ;
  int rows = rmatrix.rows() ;
  Vector currVec(rows) ;
  for (int i = 0; i < rows; ++ i) {
    currVec(i) = GetDouble( (RNumber)rmatrix.at(i, 0) ) ;
  } 
  return currVec ;
}

// get rational matrix from integer matrix (stored in double)
RMatrix Tool::GetRatioalMatrix(const Matrix& matrix) {
  int rowNum = matrix.rows() ;
  int colNum = matrix.cols() ;
  RMatrix rmatrix(rowNum, colNum) ;
  int val ;
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum; ++ j) {
      val = round( matrix(i, j) ) ;
      rmatrix.at(i, j) = val ;
    }
  }

  return rmatrix ;
}

RMatrix Tool::GetBlock(const RMatrix& matrix, 
    int startRow, int startCol, int rowNum, int colNum) {
  RMatrix newMatrix(rowNum, colNum) ;
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum; ++ j) {
      newMatrix.at(i, j) = matrix.at(startRow+i, startCol+j) ;
    }
  }
  return newMatrix ;
}

void Tool::CopyMatrix(const RMatrix& from, RMatrix& to) {
  int rowNum = from.rows() ;
  int colNum = from.cols() ;
  if (to.rows() != rowNum || to.cols() != colNum) {
    std::cerr << "Cannot copy rational matrix with different size." << std::endl ;
    std::terminate() ;
  }
  for (int i = 0; i < rowNum; ++ i) {
    for (int j = 0; j < colNum; ++ j) {
      to.at(i, j) = std::move( from.at(i, j) ) ;
    }
  }
}


void Tool::CopyMatrix(const RMatrix& from, RMatrix& to, int fromBeginRow,
      int fromBeginCol, int toBeginRow, int toBeginCol, int rowSize, int colSize) {
  if ( fromBeginRow+rowSize > from.rows() || fromBeginCol+colSize > from.cols()
      || toBeginRow+rowSize > to.rows() || toBeginCol+colSize > to.cols() ) {
    std::cerr << "Cannot copy matrix. Index is out of range." << std::endl ;
    std::terminate() ;
  }
  for (int i = 0; i < rowSize; ++ i) {
    for (int j = 0; j < colSize; ++ j) {
      to.at(i+toBeginRow, j+toBeginCol) =
          std::move( from.at(i+fromBeginRow, j+fromBeginCol) ) ;
    }
  }
}

bool Tool::IsZeroRow(const RMatrix& matrix, int rowIdx) {
  for (int j = 0; j < matrix.cols(); ++ j) {
    if ( ! matrix.at(rowIdx, j).is_zero() ) {
      return false ;
    }
  }
  return true ;
}

bool Tool::IsZeroCol(const RMatrix& matrix, int colIdx) {
  for (int i = 0; i < matrix.rows(); ++ i) {
    if ( ! matrix.at(i, colIdx).is_zero() ) {
      return false ;
    }
  }
  return true ;
}

bool Tool::IsZeroBlock(const RMatrix& matrix,
    int beginRow, int rowSize, int beginCol, int colSize) {
  for (int i = 0; i < rowSize; ++ i) {
    for (int j = 0; j < colSize; ++ j) {
      if ( ! matrix.at(beginRow+i, beginCol+j).is_zero() ) {
        return false ;
      }
    }
  }
  return true ;
}

bool Tool::ColumnsEq(const RMatrix& matrix1, const RMatrix& matrix2,
    int idx1, int idx2, bool inverse) {
  bool res = true ;
  if ( matrix1.rows() != matrix2.rows() ) {
    return false ;
  }
  // find the first coeff which is not zero
  int ratioIdx = -1 ;
  for (int i = 0; i < matrix2.rows(); ++ i) {
    if ( ! matrix2.at(i, idx2).is_zero() ) {
      ratioIdx = i ;
      break ;
    }
  }
  if (ratioIdx == -1) {
    std::cerr << "ColumnsEq(): error: zero vector." << std::endl ;
    std::terminate() ; 
  }
  RNumber zero ;
  RNumber ratio( matrix1.at(ratioIdx, idx1) / matrix2.at(ratioIdx, idx2) ) ;
  if ( (inverse && ratio > zero) || ( ! inverse && ratio < zero) ) {
    return false ;
  }
  for (int i = 0; i < matrix1.rows(); ++ i) {
    if ( matrix1.at(i, idx1) != matrix2.at(i, idx2) * ratio ) {
      res = false ;
      break ;
    }
  }

  return res ;
}

bool Tool::ColumnsEq(const RMatrix& matrix, int idx1, int idx2, bool inverse) {
  return ColumnsEq(matrix, matrix, idx1, idx2, inverse) ; 
}

ReducedMatrix Tool::GaussianElimination(const RMatrix& matrix) {
  RMatrix constraints(matrix) ;
  const int initVariNum = matrix.cols() - 1 ;
  const int constraintsCols = matrix.cols() ;
  std::vector<int> basic, dupRow ;
  // choose basic variables
  int basicIdx ;
  RNumber den, curr, ratio ;
  ReducedMatrix res( matrix.rows(), matrix.cols() ) ;

  int init_cons_num = constraints.rows() ;
  for (int i = 0; i < init_cons_num; ++ i) {
    basicIdx = -1 ;
    for (int j = 0; j < initVariNum; ++ j) {
      if ( ! constraints.at(i, j).is_zero() ) {
        basicIdx = j ;
        break ;
      }
    }
    // if basicIdx == -1, then this constraint is a duplicated one, or
    // the constraint is infeasible
    if (basicIdx != -1) {
      basic.push_back(basicIdx) ;
      res.rowIdx.push_back(i) ;
      
      den = constraints.at(i, basicIdx) ;
      for (int j = 0; j < constraintsCols; ++ j) {
        curr = constraints.at(i, j) ;
        constraints.at(i, j) = curr/den ; 
      }
      
      for (int k = 0; k < init_cons_num; ++ k) {
        if (k == i) continue ;
        ratio = constraints.at(k, basicIdx) ;
        for (int j = 0; j < constraintsCols; ++ j) {
          curr = constraints.at(k, j) ;
          constraints.at(k, j) = curr - ratio * constraints.at(i, j) ;
        } 
      }
    }
    else {
      dupRow.push_back(i) ;
    }

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  if (basicIdx == -1) {
    std::cout << "row " << i << " cannot choose basic variable." << std::endl ;
  }
  log_mtx.unlock() ;
#endif

  } 

  // infeasible
  if ( ! dupRow.empty() ) {
    for (int i : dupRow) {
      // if the constraint is infeasible
      if ( ! constraints.at(i, constraintsCols-1).is_zero() ) {
        res.succeed = false ;
        return res ;
      }
    }
  }

  for (int j = 0; j < initVariNum; ++ j) {
    if ( std::find(basic.begin(), basic.end(), j) == basic.end() ) {
      res.colIdx.push_back(j) ;
    }
  }
  for (unsigned j = 0; j < basic.size(); ++ j) {
    res.colIdx.push_back( basic[j] ) ;
  }
  CopyMatrix(constraints, res.matrix) ;
  res.basicIdx = std::move(basic) ;

#ifdef DEBUGINFO_PLP
  log_mtx.lock() ;
  std::cout << "order: " ;
  for (unsigned j = 0; j < res.colIdx.size(); ++ j) {
    std::cout << " " << res.colIdx[j] << " " ;
  }
  std::cout << std::endl ;
  std::cout << "constraints after choosing basis: "; print(res.matrix) ;
  log_mtx.unlock() ;
#endif

  return res ;
}

ReducedMatrixFloat Tool::GaussianEliminationFloat(const Matrix& matrix) {
  Matrix constraints = matrix ;
  const int initVariNum = matrix.cols() - 1 ;
  const int constraintsCols = matrix.cols() ;
  std::vector<int> basic, dupRow ;
  // choose basic variables
  int basicIdx ;
  double den, curr, ratio ;
  ReducedMatrixFloat res ;
 
  int init_cons_num = constraints.rows() ;
  for (int i = 0; i < init_cons_num; ++ i) {
    basicIdx = -1 ;
    for (int j = 0; j < initVariNum; ++ j) {
      if ( constraints(i,j) != 0 ) {
        basicIdx = j ;
        break ;
      }
    }
    // if basicIdx == -1, then this constraint is a duplicated one, or
    // the constraint is infeasible
    if (basicIdx != -1) {
      basic.push_back(basicIdx) ;
      res.rowIdx.push_back(i) ;
      
      den = constraints(i, basicIdx) ;
      for (int j = 0; j < constraintsCols; ++ j) {
        curr = constraints(i, j) ;
        constraints(i, j) = curr/den ; 
      }
      
      for (int k = 0; k < init_cons_num; ++ k) {
        if (k == i) continue ;
        ratio = constraints(k, basicIdx) ;
        for (int j = 0; j < constraintsCols; ++ j) {
          curr = constraints(k, j) ;
          constraints(k, j) = curr - ratio * constraints(i, j) ;
          // if the value is small, set as zero for eliminating dup constraint
          if (constraints(k, j) < 1e-6 && constraints(k, j) > -1e-6) {
            constraints(k, j) = 0 ;
          }
        } 
      }

      // set 0 or 1 for floating point
      constraints(i, basicIdx) = 1 ;
      for (int k = 0; k < init_cons_num; ++ k) {
        if (k == i) continue ;
        constraints(k, basicIdx) = 0 ;
      }
    }
    else {
      dupRow.push_back(i) ;
    }
  } 

  // infeasible
  if ( ! dupRow.empty() ) {
    for (int i : dupRow) {
      // if the constraint is infeasible
      if (constraints(i, constraintsCols-1) != 0) {
        res.succeed = false ;
        return res ;
      }
    }
  }

  for (int j = 0; j < initVariNum; ++ j) {
    if ( std::find(basic.begin(), basic.end(), j) == basic.end() ) {
      res.colIdx.push_back(j) ;
    }
  }
  for (unsigned j = 0; j < basic.size(); ++ j) {
    res.colIdx.push_back( basic[j] ) ;
  }
  res.matrix = constraints ;
  res.basicIdx = std::move(basic) ;

  return res ;
}


void Tool::BinomialCoeff(const std::vector<int>& idx, unsigned k, std::vector< std::vector<int> >& res, std::list<int> prev) { 
  unsigned num = idx.size() ;
  if (k == 0) {
    std::vector<int> curr ;
    for (auto it = prev.begin(); it != prev.end(); ++ it) {
      curr.push_back(*it) ;
    }
    std::sort( curr.begin(), curr.end() ) ;
    res.push_back( std::move(curr) ) ;
    return ;
  }
  if (k == num) {
    std::vector<int> curr ;
    for (auto it = prev.begin(); it != prev.end(); ++ it) {
      curr.push_back(*it) ;
    }

    for (int i : idx) {
      curr.push_back(i) ;
    }
    std::sort( curr.begin(), curr.end() ) ;
    res.push_back( std::move(curr) ) ;
    return ;
  }
  for (unsigned i = 0; i < num; ++ i) {
    if (num - i < k) break ;
    prev.push_front(idx[i]) ;
    std::vector<int> idx2 ;
    for (unsigned j = i+1; j < num; ++ j) {
      idx2.push_back(idx[j]) ;
    }
    BinomialCoeff(idx2, k-1, res, prev) ;
    prev.pop_front() ;
  }
}


std::vector<int> Tool::GetNonDupIdx(const Matrix& m, bool useFloat) {
  bool dup ;
  std::vector<int> idx ;
  if ( ! m.row(0).isZero() ) {
    idx.push_back(0) ;
  }
  int consNum = m.rows() ;
  int colNum = m.cols() ;
  for (int i = 1; i < consNum; ++ i) {
    dup = false ;
    if ( m.row(i).isZero() ) {
      continue ;
    }
    for (int k = i-1; k >= 0; -- k) {
      if (useFloat) { 
        // ratio will must be assigned later
        double curr, ratio = 0 ;
        bool allEq = true ;
        for (int j = 0; j < colNum-1; ++ j) {
          if (m(k,j) != 0) {
            ratio = m(i,j) / m(k,j) ;
          }
        }
        for (int j = 0; j < colNum; ++ j) {
          if (m(k,j) == 0) {
            if (m(i,j) == 0) {
              continue ;
            }
            else {
              allEq = false ;
              break ;
            }
          }
          curr = m(i,j) / m(k,j) ;
          // ratio must have been initialized.
          // Otherwise either m.row(k) and m.row(i) are all zeros, which is not
          // possible, or allEq == flase and break.
          if ( curr != ratio ) {
            allEq = false ;
            break ;
          }
        } 
        if (allEq) {
          dup = true ;
          break ;
        }
      }
      // the coeffs have been divided by gcd
      else {
        if ( m.row(i) == m.row(k) ) {
          dup = true ;
          break ;
        }
      }
    }
    if (! dup) {
      idx.push_back(i) ;
    }
  }
  return idx ;
}

std::vector<int> Tool::GetNonDupIdx(const RMatrix& m) {
  bool dup ;
  std::vector<int> idx ;
  if ( ! Tool::IsZeroRow(m,0) ) {
    idx.push_back(0) ;
  }
  int consNum = m.rows() ;
  int colNum = m.cols() ;
  RNumber curr, ratio ;
  for (int i = 1; i < consNum; ++ i) {
    dup = false ;
    if ( Tool::IsZeroRow(m, i) ) {
      continue ;
    }
    for (int k = i-1; k >= 0; -- k) {
      bool allEq = true ;
      for (int j = 0; j < colNum-1; ++ j) {
        if ( ! m.at(k,j).is_zero() ) {
          ratio = m.at(i,j) / m.at(k,j) ;
        }
      }
      for (int j = 0; j < colNum; ++ j) {
        if ( m.at(k,j).is_zero() ) {
          if ( m.at(i,j).is_zero() ) {
            continue ;
          }
          else {
            allEq = false ;
            break ;
          }
        }
        curr = m.at(i,j) / m.at(k,j) ;
        if ( curr != ratio ) {
          allEq = false ;
          break ;
        }
      } 
      if (allEq) {
        dup = true ;
        break ;
      }
    }
    if (! dup) {
      idx.push_back(i) ;
    }
  }
  return idx ;
}

/*******************************tmp code********************************/
bool Tool::CmpFromList(const ChooseFrom& vec1, const ChooseFrom& vec2) { 
  assert( vec1.currIdx.size() == vec2.currIdx.size() ) ;
  std::vector<int> v1(vec1.currIdx) ;
  std::vector<int> v2(vec2.currIdx) ;
  std::sort( v1.begin(), v1.end() ) ;
  std::sort( v2.begin(), v2.end() ) ;
  for (unsigned i = 0; i < v1.size(); ++ i) {
    if ( v1[i] > v2[i] ) return false ;
    if ( v1[i] < v2[i] ) return true ;
  }
  return true ;
}


bool Tool::LexicoLess(const RMatrix& vec1, const RMatrix& vec2) {
  assert( vec1.cols() == vec2.cols() ) ;
  for (unsigned j = 0; j < vec1.cols(); ++ j) {
    if ( vec1.at(0,j) < vec2.at(0,j) ) return true ;
    if ( vec1.at(0,j) > vec2.at(0,j) ) return false ;
  }
  return true ;
}


bool Tool::LexicoPositive(const RMatrix& matrix, int rowIdx) {
  RNumber zero ;
  for (unsigned j = 0; j < matrix.cols(); ++ j) {
    if (matrix.at(rowIdx, j) > zero) return true ;
    if (matrix.at(rowIdx, j) < zero) return false ;
  }
  return false ;
}

/***********************************************************************/
