#include <iostream>
#include "pplp_apronInterface.h"
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"
#include "pplp_tool.h"
#ifdef _TBB
#include "tbb/task_scheduler_init.h"
#elif defined(_OPENMP)
#include <omp.h>
#endif

using namespace PPLP ;

int main (int argc, char* argv[]) {
  if(argc != 2 && argc != 3) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <filename> [nb_threads]" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ;
  bool res, allequal = true ;
  int idx = 0 ;
  std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;

#ifdef _TBB
  tbb::task_scheduler_init scheduler_init(argc>=3 ? atoi(argv[2]) : 1);
#endif
#ifdef _OPENMP
  omp_set_dynamic(0) ; 
  omp_set_num_threads(argc>=3 ? atoi(argv[2]) : 1) ;
#endif

  int curr ;
  for (auto& poly : polyVec) {
    std::cout << polyName << ":" << idx << std::endl ;
    if (poly.get_eq_constraint_num() != 0) {
      std::cout << "Cannot compare Polyhedron having equalities with apron" << std::endl ;
      continue ;
    }

    poly.Minimize() ; 
    // all the coeffs are integers, convert it into rational
    // and then compare with Apron
    std::vector<int> activeIdx = poly.GetActiveIdx() ;
    int consNum = poly.get_constraint_num() ;
    int subConsNum = poly.GetActiveIdx().size() ;
    int variNum = poly.get_variable_num() ; 
    RMatrix rationalPoly(subConsNum, variNum + 1) ;
    for (int i = 0; i < subConsNum; ++ i) {
      int currIdx = activeIdx[i] ;
      for (int j = 0; j < variNum; ++ j) { 
        curr = - poly.GetCoef(currIdx, j) ;
        rationalPoly.at(i, j) = curr ;
      }
      curr = poly.GetConstant(currIdx) ;
      rationalPoly.at(i, variNum) = curr ;
    }
    ApronInterface raytracingMini(rationalPoly) ;

    RMatrix rationalApron(consNum, variNum + 1) ;
    for (int i = 0; i < consNum; ++ i) {
      for (int j = 0; j < variNum; ++ j) { 
        curr = - poly.GetCoef(i, j) ;
        rationalApron.at(i, j) = curr ;
      }
      curr = poly.GetConstant(i) ;
      rationalApron.at(i, variNum) = curr ;
    }
    ApronInterface apronMini(rationalApron) ; 

    if (apronMini.GetConsNum() != subConsNum) {
      res = false ; 
    }
    else {
      res = raytracingMini.CompareAbstract(apronMini) ;
    }

    if (res == true) {
      //std::cout << "Minimized polyhedra " << ioInter.n << ": " 
      //      << idx << " are equal." << std::endl ;
    }
    else {
      std::cout << "Minimized polyhedra " << polyName << ": " 
          << idx << " are NOT equal." << std::endl ;
      std::cout << "raytracing:" << std::endl ;
      raytracingMini.Print() ;
      std::cout << "Apron:" << std::endl ;
      apronMini.Print() ;
    }

      //allequal = false ;
      //std::terminate() ;
    ++ idx ;
  } 
  /*
  if (allequal) {
    std::cout << "All the minimized polyhedra are equal." << std::endl ;    
  }
  */
  return 0 ;
}
