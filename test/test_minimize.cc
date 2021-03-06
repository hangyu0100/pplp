#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"
#ifdef _TBB
#include "tbb/task_scheduler_init.h"
#elif defined(_OPENMP)
#include <omp.h>
#endif

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc < 2 || argc > 4) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <filename> [nb_threads] [output path]" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter ;
  std::ofstream ofs ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
  if (argc == 4) {
    std::string outpath = std::string(argv[3]) ;
    ofs.open(outpath, std::fstream::out | std::fstream::app) ; 
    if ( ! ofs.is_open() ) {
      std::cerr << "Cannot open file." << std::endl ;
      return -1 ;
    }
  }

#ifdef _TBB
  tbb::task_scheduler_init scheduler_init(argc>=3 ? atoi(argv[2]) : 1);
#endif
#ifdef _OPENMP
  omp_set_dynamic(0) ; 
  omp_set_num_threads(argc>=3 ? atoi(argv[2]) : 1) ;
#endif

  //double totalTime = 0 ;
  int idx = 0 ;
  //std::cout << "P_" << ioInter.get_cons_num() << "_" << ioInter.get_redundancy()
  //    << "_" << ioInter.get_vari_num() << "_" << ioInter.get_zero_num() 
  //    << ".poly" << std::endl ; 
  
  std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;

  for (auto& poly : polyVec) {
    std::cout << polyName << ": " << idx << std::endl ;
    // polyNum, total consNum, redundancy, variNum, zeros, generators, id 
    auto start = std::chrono::steady_clock::now() ;
    Polyhedron minimizedP = poly.GetMinimizedPoly(true) ;
    auto end = std::chrono::steady_clock::now() ;
    std::chrono::duration<double> diff = (end - start) * 1000 ;
    //std::cout << "Running time: " << diff.count() << " ms" << std::endl ;
   
    minimizedP.Print() ;
   
  }
  
  return 0 ;
} 
