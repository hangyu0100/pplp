#include <cstdlib>
#include <iostream>
#include <fstream>
#include "pplp_ioInterface.h"
#include "pplp_plp.h"
#include "pplp_tbbParallel.h"
#include "pplp_apronInterface.h"

#ifdef VERIMAG_POLYHEDRA_PLP_TBB
#include "tbb/task_scheduler_init.h"
#endif

#include <chrono>

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc < 2 || argc > 4) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <filename> [nb_threads] [nb_initial_points]" << std::endl ;
    return -1 ; 
  }

#ifdef VERIMAG_POLYHEDRA_PLP_TBB
  tbb::task_scheduler_init scheduler_init( argc>=3 ? atoi(argv[2]) : 1);
#endif
  int nb_initial_points = argc>=4 ? atoi(argv[3]) : 1;
  
  //open the log file
  //std::string outpath = "./plp_time.log" ;
  /*
  std::ofstream ofs(outpath, std::fstream::out | std::fstream::app) ; 
  if ( ! ofs.is_open() ) {
    std::cerr << "Cannot open file." << std::endl ;
    return -1 ;
  }
  */
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
  std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;

  double total = 0 ;
  bool res = true ;
  int idx = 0, curr ;
  for (auto& poly : polyVec) {
    std::cout << polyName << ":" << idx << std::endl ;
    //auto start = std::chrono::steady_clock::now() ;
    TbbParallel tbb(poly) ;
    tbb.PlpParallel(nb_initial_points) ;
    //auto end = std::chrono::steady_clock::now() ;
    //std::chrono::duration<double> diff = (end - start) * 1000 ;
    //total += diff.count() ;
    //Polyhedron plpMiniPoly = tbb.GetOptPoly() ;

    //plpMiniPoly.Print() ;

    RMatrix tbbOpt = tbb.GetOptimalMatrix() ;
    ApronInterface plpMini(tbbOpt) ;

    int consNum = poly.get_constraint_num() ;
    int variNum = poly.get_variable_num() ; 
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

    int plpConsNum = tbbOpt.rows() ;
    if (apronMini.GetConsNum() != plpConsNum) {
      res = false ; 
    }
    else {
      res = plpMini.CompareAbstract(apronMini) ;
    }

    if (res == true) {
      //std::cout << "Minimized polyhedra " << ioInter.n << ": " 
      //      << idx << " are equal." << std::endl ;
    }
    else {
      std::cout << "Minimized polyhedra " << polyName << ": " 
          << idx << " are NOT equal." << std::endl ;
      std::cout << "plp:" << std::endl ;
      plpMini.Print() ;
      std::cout << "Apron:" << std::endl ;
      apronMini.Print() ;
    }

      //allequal = false ;
      //std::terminate() ;
    ++ idx ;
  } 

  //std::cout << "total time: " << total << std::endl ;
  //std::cout << "spin count: " << spin_count << std::endl;
  return 0 ;
} 
