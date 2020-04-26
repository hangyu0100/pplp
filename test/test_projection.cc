#include "pplp_plp.h"
#include "pplp_ioInterface.h"
#include "pplp_tbbParallel.h"
#ifdef _TBB
#include "tbb/task_scheduler_init.h"
#elif defined(_OPENMP)
#include <omp.h>
#endif
#include "pplp_apronInterface.h"
//#include <flint.h>
#include <flint/flint.h>

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc < 2 || argc > 5) {
    std::cerr << "usage: " << argv[0] 
        << " <filename> [number of projected variables]"
        << " [nb_threads] [nb_initial_points]" << std::endl ;
    return -1 ; 
  }

  {
    IoInterface ioInter ;
    std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
#ifdef NO_STRING
    std::string polyName = "TOTO";
#else
    std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;
#endif

    int projNum = argc >= 3 ? atoi(argv[2]) : (ioInter.get_vari_num() / 2) ; 
    if ( projNum >= ioInter.get_vari_num() ) {
      std::cerr << "Cannot project " << projNum << " variables" << std::endl ;
      std::terminate() ;
    }
    
#ifdef _TBB 
    tbb::task_scheduler_init scheduler_init( argc>=4 ? atoi(argv[3]) : 1);
#endif
#ifdef _OPENMP
  omp_set_dynamic(0) ; 
  omp_set_num_threads(argc>=4 ? atoi(argv[3]) : 1) ;
#endif
    
  int nb_initial_points = argc>=5 ? atoi(argv[4]) : 1;
  
  double total = 0.0;
  int id = 0 ;
  for (auto& poly : polyVec) {
    //std::cout << polyName << "_" << id << std::endl ;
    //Polyhedron miniPoly = poly.GetMinimizedPoly() ;
    Polyhedron miniPoly = poly.GetMinimizedPoly(false) ;
    auto start = std::chrono::steady_clock::now() ;
    TbbParallel tbb(miniPoly) ;
    tbb.PlpParallel(projNum, nb_initial_points) ;
    auto end = std::chrono::steady_clock::now() ;
    std::chrono::duration<double> diff = (end - start) * 1000 ;

    std::cout << diff.count() << std::endl ;
      
    //total += diff.count() ;
    ++ id ;
  }
  //std::cout << "total time: " << total << " work time:" << work_microseconds*1E-3 << " " << "spin count: " << spin_count << std::endl;
  }
#ifdef VERIMAG_POLYHEDRA_PLP_OPENMP
#pragma omp parallel
#endif
  {
    glp_free_env();
    flint_cleanup();
  }
  return 0 ;
}
