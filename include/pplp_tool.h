#ifndef _RAYTRACING_TOOL
#define _RAYTRACING_TOOL

#include <list>
#include <eigen3/Eigen/Dense>
#include <flint/fmpq_matxx.h>
typedef Eigen::MatrixXd Matrix ;
typedef Eigen::RowVectorXd Vector ;
typedef Eigen::MatrixXi MatrixZ ;
typedef Eigen::RowVectorXi VectorZ ;
typedef flint::fmpqxx RNumber ;
typedef flint::fmpq_matxx RMatrix ;
typedef flint::fmpzxx ZNumber ;
#include "pplp_double.h"

namespace PPLP {

// for Gaussian Elimination
struct ReducedMatrix {
  ReducedMatrix(int rowNum, int colNum) : matrix(rowNum, colNum), succeed(true) {} ;
  RMatrix matrix ;
  std::vector<int> rowIdx ;
  std::vector<int> colIdx ;
  std::vector<int> basicIdx ;
  bool succeed ;
} ;

struct ReducedMatrixFloat {
  ReducedMatrixFloat() : succeed(true) {} ;
  Matrix matrix ;
  std::vector<int> rowIdx ;
  std::vector<int> colIdx ;
  std::vector<int> basicIdx ;
  bool succeed ;
} ;

enum RegionState {zero=-4, noCons=-3, flat=-2, normal=-1} ;

/*******************************tmp code********************************/

struct ConstraintIdx {
    ConstraintIdx() : regIdx(-1), frontierIdx(-1) {}
    ConstraintIdx(int r, int f) : regIdx(r), frontierIdx(f) {} 
    int regIdx ;
    int frontierIdx ;
  } ;
  
struct ChooseFrom {
    ChooseFrom(const std::vector<int>& curr, const RMatrix matrix,
        const ConstraintIdx& f) : currIdx(curr), pertMatrix(matrix), from(f) {}
    std::vector<int> currIdx ;
    RMatrix pertMatrix ;
    ConstraintIdx from ;
  } ;

/*******************************tmp code ends***************************/






class Tool {
public: 
  static long int GetGcd(const long int val1, const long int val2) {
    long int res = val2 == 0 ? val1 : GetGcd(val2, val1 % val2) ;
    return std::abs(res) ;
  } 

  static double GetGcd(const double val1, const double val2) {
    double res = val2 == 0 ? val1 : GetGcd( val2, fmod(val1, val2) ) ;
    return std::abs(res) ;
  }

  static double GetGcd(const Vector& cons, double constant) {
    double gcd = Tool::GetGcd( constant, cons(0) ) ;
    if (gcd != 1) {
      for (int j = 1; j < cons.cols(); ++ j) {
        gcd = Tool::GetGcd( gcd, cons(j) ) ;
        if (gcd == 1) break ;
      }
    }
    return gcd ;
  }

  static double GetDotProductThreshold(const Vector& vec1, const Vector& vec2) ;
  //static double GetDotProductThreshold(int size, double norm1, double norm2) ;

  static double GetDouble(RNumber rational) {
    double num = rational.num().to<double>() ;
    double den = rational.den().to<double>() ;
    return num / den ;
  }
  static Matrix GetDoubleMatrix(const RMatrix& rmatrix) ;
  static Matrix GetDoubleBlock(const RMatrix& rmatrix, int beginRow, int beginCol,
      int rowNum, int colNum) ;
  static Vector GetDoubleVector(const RMatrix& rmatrix) ;

  static RMatrix GetRatioalMatrix(const Matrix& matrix) ;

  static RMatrix GetBlock(const RMatrix& matrix, 
      int startRow, int startCol, int rowNum, int colNum) ;

  static void CopyMatrix(const RMatrix& from, RMatrix& to) ;
  static void CopyMatrix(const RMatrix& from, RMatrix& to, int fromBeginRow,
      int fromBeginCol, int toBeginRow, int toBeginCol, int rowSize, int colSize) ;

  static int GetIntOfFloat(double val) {
    int valInt = val > 0 ? val + 0.5 : val - 0.5 ;
    return valInt ;
  }

  static bool IsZeroRow(const RMatrix& matrix, int rowIdx) ;
  static bool IsZeroCol(const RMatrix& matrix, int colIdx) ;
  static bool IsZeroBlock(const RMatrix& matrix,
      int beginRow, int rowSize, int beginCol, int colSize) ;

  static bool ColumnsEq(const RMatrix& matrix1, const RMatrix& matrix2,
    int idx1, int idx2, bool inverse=false) ;

  static bool ColumnsEq(const RMatrix& matrix,
    int idx1, int idx2, bool inverse=false) ;

  // Gaussian Elimination
  static ReducedMatrix GaussianElimination(const RMatrix& matrix) ;
  static ReducedMatrixFloat GaussianEliminationFloat(const Matrix& matrix) ;

  static void BinomialCoeff(const std::vector<int>& idx, unsigned k,
      std::vector< std::vector<int> >& res, std::list<int> prev = std::list<int>()) ;

  static std::vector<int> GetNonDupIdx(const Matrix& m, bool useFloat = false) ;
  static std::vector<int> GetNonDupIdx(const RMatrix& m) ;


/*******************************tmp code********************************/
  static bool CmpFromList(const ChooseFrom& vec1, const ChooseFrom& vec2) ; 

  static bool LexicoLess(const RMatrix& vec1, const RMatrix& vec2) ;

  static bool LexicoPositive(const RMatrix& matrix, int rowIdx) ;

/***********************************************************************/

  

} ;

}

#endif
