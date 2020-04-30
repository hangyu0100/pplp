/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
 * GlpkInterface class provides functions to access GLPK library.
 * _glp is used to store the result of simplex, and the result can be used in 
 * the next process of solving simplex. Currently it is just used 
 * in GetSatPoint(), as GetSatPoint() add the constraints incrementally.
*******************************************************************************/

#ifndef _RAYTRACING_GLPKINTER
#define _RAYTRACING_GLPKINTER

#include <glpk.h>
#include <eigen3/Eigen/Dense>
#include "pplp_polyhedron.h"


#if defined(DEBUGINFO_RAYTRACING) || defined(PRINT_WARNING)
#include <mutex>
extern std::mutex log_mtx_raytracing;
#endif

namespace PPLP {
class GlpkInterface {
public:
  GlpkInterface() ;
  ~GlpkInterface() ;
  static Point GetCentralPoint(Polyhedron& poly, bool checkFlat=false) ;
  static Point GetSatPoint(Polyhedron& poly, double epsilon = _epsilon) ;
  Point GetIrddWitness(const std::vector<int>& headIdx, const Polyhedron& poly,
      bool cone, double epsilon = _epsilon) ;
  static bool Sat(const Polyhedron& poly, int idx, double threshold = _threshold) ;
  bool Simplex( const Polyhedron& poly, const Vector& obj, 
      bool variNonNeg = false, bool askFeasible = false,
      bool getBasis = true, std::vector<int> nonNegIdx = std::vector<int>() ) ;
  double get_simplex_optimal() ;
  void GetBasis() ;
  const std::vector<int>& get_basic_idx() ; 
  const std::vector<int>& get_non_basic_idx() ;
  int GetVariState(int idx) ;
  double GetVariVal(int idx) ;
  int GetRowNum() ;
  int GetColNum() ;
  static double get_threshold() {
    return _threshold ;
  }
  void PrintValue() ;
  bool is_result_unsure() {
    return _result_unsure ;
  }
private:
  static double _epsilon ;  
  static double _threshold ;  
  glp_prob* _glp ;
  std::vector<int> _basic_idx ;
  std::vector<int> _non_basic_idx ;
  double _simplex_optimal ;
  Vector _first_constraint ;
  int _currIdx ;
  bool _result_unsure ;
} ;

}

#endif
