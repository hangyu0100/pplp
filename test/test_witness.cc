#include <iostream>
#include <vector>
#include "pplp_ioInterface.h"

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <filename>" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
  int currIdx ;
  for (auto& poly : polyVec) {
    poly.Minimize() ;
    std::cout << "the internal point is:" << std::endl 
      << poly.get_internal_point() << std::endl ;
    const std::vector<Point>& witness = poly.GetWitness() ;
    for (int i = 0; i < (int)poly.GetActiveIdx().size(); ++ i) {
      currIdx = poly.GetActiveIdx().at(i) ;
      std::cout << "constraint " << currIdx << " is:" << std::endl ;
      poly.PrintConstraint(currIdx) ;
      std::cout << "witness point is:" << std::endl 
        << witness[i].get_coordinates() << std::endl ;
    } 
  }
  return 0 ;
} 
