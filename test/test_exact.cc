#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"

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
  for (auto& poly : polyVec) { 
    //poly.Print() ;
    
    //VectorZ obj(2) ;
    //obj[0] = -1 ;
    //obj[1] = -1 ;
    
    bool res = poly.GetExactSolution() ;
    //bool res = poly.GetExactSolution(obj) ;
    if (res == false) {
      continue ;
    } 

    std::cout << std::endl << std::endl ;
  }

  return 0 ;
} 
