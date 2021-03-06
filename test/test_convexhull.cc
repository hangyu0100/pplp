#include "pplp_plp.h"
#include "pplp_ioInterface.h"
#include "pplp_tbbParallel.h"
#ifdef _TBB
#include "tbb/task_scheduler_init.h"
#elif defined(_OPENMP)
#include <omp.h>
#endif

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc < 3 || argc > 5) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename 
        << " <filename1> <filename2>"
        << " [nb_threads] [nb_initial_points]" << std::endl ;
    return -1 ; 
  }

  IoInterface ioInter1 ;
  std::vector<Polyhedron> polyVec1 = ioInter1.LoadPolyhedra(argv[1]) ; 
  IoInterface ioInter2 ;
  std::vector<Polyhedron> polyVec2 = ioInter2.LoadPolyhedra(argv[2]) ; 
  std::string polyName1 = "P_" + std::to_string( ioInter1.get_cons_num() )
      + "_" + std::to_string( ioInter1.get_redundancy() )
      + "_" + std::to_string( ioInter1.get_vari_num() )
      + "_" + std::to_string( ioInter1.get_zero_num() ) ;
  std::string polyName2 = "P_" + std::to_string( ioInter2.get_cons_num() )
      + "_" + std::to_string( ioInter2.get_redundancy() )
      + "_" + std::to_string( ioInter2.get_vari_num() )
      + "_" + std::to_string( ioInter2.get_zero_num() ) ;

#ifdef VERIMAG_POLYHEDRA_PLP_TBB
  tbb::task_scheduler_init scheduler_init( argc>=4 ? atoi(argv[3]) : 1);
#endif
#ifdef _OPENMP
  omp_set_dynamic(0) ; 
  omp_set_num_threads(argc>=4 ? atoi(argv[3]) : 1) ;
#endif


  int nb_initial_points = argc>=5 ? atoi(argv[4]) : 1;
  
  int polyNum = polyVec1.size() < polyVec2.size() ?
      polyVec1.size() : polyVec2.size() ;

  for (int i = 0; i < polyNum; ++ i) {
    //std::cout << polyName1 << ", " << polyName2 << " : " << i << std::endl ;
    Polyhedron miniPoly1 = polyVec1[i].GetMinimizedPoly() ;
    Polyhedron miniPoly2 = polyVec2[i].GetMinimizedPoly() ;
    TbbParallel tbb(miniPoly1, miniPoly2) ;
    auto start_plp = std::chrono::high_resolution_clock::now();
    tbb.PlpConvexHull(nb_initial_points) ;
    auto end_plp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_plp = (end_plp - start_plp) * 1e3 ;

    Polyhedron projPoly = tbb.GetOptPoly() ;
    projPoly.Print() ;

  }
  return 0 ;
}
