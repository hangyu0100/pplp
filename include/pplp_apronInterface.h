/*******************************************************************************
 * Copyright (c) 2017 Jan. Verimag. All rights reserved.
 * @author Hang YU
 * This class provides some functions used library APRON.
 * An object of ApronInterface stores the constraints of a polyhedron. 
 * When we transform ap_lincons1_array_t to ap_abstract1_t, it will eliminate 
 * the redundant constraints, so we do not need to do extra work for 
 * minimization.
*******************************************************************************/

#ifndef _RAYTRACING_APRONINTER
#define _RAYTRACING_APRONINTER

#include <unordered_set>
#include <ap_global1.h>
#include <pk.h>
#include "pplp_polyhedron.h"

namespace PPLP {

class ApronInterface {
public:
  ApronInterface(const Polyhedron& poly) ;
  ApronInterface(const RMatrix& optimals) ;
  ~ApronInterface() ;
  Polyhedron GetPoly() ;
  Polyhedron GetPoly (ap_abstract1_t abs) ;
  void GetAbstractVal() ;
  bool CompareAbstract(const ApronInterface& apron) ;
  int GetGeneratorNum() ;
  ap_abstract1_t get_abs() const ;
  void Print() ;
  // TODO get the list of vetices
  std::vector<double> get_gen_array() ;
  void PrintGenerator() ;
  double GetGenVal(int idx, int dim) ;
  bool CmpProjectedPoly(const std::vector<int>& idx, const RMatrix& matrix) ;
  void Project(int projNum) ;
  bool CmpConvexHull(const Polyhedron& poly2, const RMatrix& matrix) ;
  ap_lincons1_array_t ConvexHull(const Polyhedron& poly2) ;
  int GetConsNum() ;
private:  
  int _vari_num ;
  char** _name_of_vari ;
  ap_manager_t* _man ; 
  ap_environment_t* _env ;
  ap_lincons1_array_t _cons_array ;
  ap_abstract1_t _abs ;
  ap_generator1_array_t _generator ;
  // store the generators as a 1D array, the length is dimension*size
  std::vector<double> _gen_array ;
  void FlintToApron(ap_scalar_t* scalar, const RNumber& val) ;
} ;

}

#endif
