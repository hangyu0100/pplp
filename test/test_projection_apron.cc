#include "pplp_plp.h"
#include "pplp_ioInterface.h"
#include "pplp_tbbParallel.h"
#ifdef _TBB
#include "tbb/task_scheduler_init.h"
#elif defined(_OPENMP)
#include <omp.h>
#endif
#include "pplp_apronInterface.h"

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc < 2 || argc > 5) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename 
        << " <filename> [number of projected variables]"
        << " [nb_threads] [nb_initial_points]" << std::endl ;
    return -1 ; 
  }
  
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
  std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;
  
  int projNum = argc >= 3 ? atoi(argv[2]) : (ioInter.get_vari_num() / 2) ; 
  if ( projNum >= ioInter.get_vari_num() ) {
    std::cerr << "Cannot project " << projNum << " variables" << std::endl ;
    std::terminate() ;
  }

#ifdef _TBB
  tbb::task_scheduler_init scheduler_init(argc>=4 ? atoi(argv[3]) : 1);
#endif
#ifdef _OPENMP
  omp_set_dynamic(0) ; 
  omp_set_num_threads(argc>=4 ? atoi(argv[3]) : 1) ;
#endif
  int nb_initial_points = argc>=5 ? atoi(argv[4]) : 1;

  for (auto& poly : polyVec) {
    std::cout << polyName << " : " << poly.get_id() << std::endl ;
    Polyhedron miniPoly = poly.GetMinimizedPoly() ;
    //miniPoly.Print() ;
    auto start_apron = std::chrono::high_resolution_clock::now();
    ApronInterface apronInter(poly) ;
    apronInter.Project(projNum) ;
    auto end_apron = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_apron = (end_apron - start_apron) * 1e3 ;
    TbbParallel tbb(miniPoly) ;
    auto start_plp = std::chrono::high_resolution_clock::now();
    tbb.PlpParallel(projNum, nb_initial_points) ;
    auto end_plp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_plp = (end_plp - start_plp) * 1e3 ;

    std::cout << "apron," << diff_apron.count() << ",plp," << diff_plp.count() << std::endl ;

    /*
    if ( apronInter.CmpProjectedPoly( tbb.get_proj_idx(), tbb.GetOptimalMatrix() ) ) {
      std::cout << "*******************Projection results in " << polyName << " of poly " 
          << poly.get_id() << " are the same." << std::endl ;
    }
    else {
      Polyhedron projPoly = tbb.GetOptPoly() ;
      std::cout << "*******************Projection results in " << polyName << " of poly " 
          << poly.get_id() << " are NOT the same." << std::endl ;
      //projPoly.Print() ;
    }
    */
  }
  return 0 ;
}
