#include <iostream>
#include <vector>
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cerr << "usage: " << argv[0] << " <inside_poly_filename> "
      << "<outside_poly_filename>" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter1, ioInter2 ;
  std::vector<Polyhedron> polyVec1 = ioInter1.LoadPolyhedra(argv[1]) ; 
  std::vector<Polyhedron> polyVec2 = ioInter2.LoadPolyhedra(argv[2]) ; 

  // test if miniPoly1 is included in miniPoly2
  int polyNum = polyVec1.size() ;
  for (int i = 0; i < polyNum; ++ i) {
    Polyhedron miniPoly1 = polyVec1[i].GetMinimizedPoly() ;
    Polyhedron miniPoly2 = polyVec2[i].GetMinimizedPoly() ;
    if (miniPoly1.IsEqualTo(std::move( miniPoly2) ) == true) {
      std::cout << "The two polyhedra are equal." << std::endl ;
    }
    else {
      std::cout << "The two polyhedra are NOT equal." << std::endl ; 
    }
  }
  return 0 ;

}
