#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc != 4) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <inside_poly_filename> "
      << "<outside_poly_filename> <filepath>" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter1, ioInter2 ;
  std::vector<Polyhedron> polyVec1 = ioInter1.LoadPolyhedra(argv[1]) ; 
  std::vector<Polyhedron> polyVec2 = ioInter2.LoadPolyhedra(argv[2]) ; 
  
  std::string filepath(argv[3]) ;
  std::ofstream ofs(filepath, std::fstream::out | std::fstream::app) ; 
  if ( ! ofs.is_open() ) {
    std::cerr << "Cannot open file." << std::endl ;
    return -1 ;
  }
  ofs << "Raytracing,Standard,Include" << std::endl ;
  
  // test if miniPoly1 is included in miniPoly2
  bool res1, res2 ;
  int polyNum = polyVec1.size() ;
  for (int i = 0; i < polyNum; ++ i) {
    Polyhedron miniPoly1 = polyVec1[i].GetMinimizedPoly() ;
    Polyhedron miniPoly2 = polyVec2[i].GetMinimizedPoly() ;
    std::cout << "Raytracing:" << std::endl ;
    auto start = std::chrono::steady_clock::now() ;
    res1 = miniPoly2.Include(miniPoly1) ;
    auto end = std::chrono::steady_clock::now() ;
    std::chrono::duration<double> diff = (end - start) * 1000 ;
    std::cout << "Running time: " << diff.count() << " ms" << std::endl ;
    ofs << diff.count() << "," ;
    if (res1 == true) {
      std::cout << "The first polyhedron is included in the second one" << std::endl ;
    }
    else {
      std::cout << "The first polyhedron is NOT included in the second one" << std::endl ; 
    }
    std::cout << "Standard:" << std::endl ;
    start = std::chrono::steady_clock::now() ;
    res2 = miniPoly2.IncludeStandard(miniPoly1) ;
    end = std::chrono::steady_clock::now() ;
    diff = (end - start) * 1000 ;
    std::cout << "Running time: " << diff.count() << " ms" << std::endl ;
    ofs << diff.count() << "," ;
    if (res2 == true) {
      std::cout << "The first polyhedron is included in the second one" << std::endl ;
    }
    else {
      std::cout << "The first polyhedron is NOT included in the second one" << std::endl ; 
    }
    std::cout << std::endl ;
    if (res1 != res2) {
      std::cerr << "The two method got different results." << std::endl ;
      std::terminate() ;
    }
    ofs << res1 << std::endl ;
  }
  return 0 ;
}
