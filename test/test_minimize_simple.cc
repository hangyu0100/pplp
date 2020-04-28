#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <chrono>
#include "pplp_ioInterface.h"
#include "pplp_glpkInterface.h"

using namespace PPLP ;

int main(int argc, char* argv[]) {
  if(argc != 3 && argc != 2) {
    std::string filepath = argv[0] ;
    std::size_t found = filepath.find_last_of("/") ;
    std::string filename = filepath.substr(found+1) ;
    std::cerr << "usage: " << filename << " <filename> [output path]" << std::endl ;
    return -1 ; 
  }
  IoInterface ioInter ;
  std::vector<Polyhedron> polyVec = ioInter.LoadPolyhedra(argv[1]) ; 
  std::ofstream ofs ;
  if (argc == 3) {
    std::string outpath = std::string(argv[2]) ;
    ofs.open(outpath, std::fstream::out | std::fstream::app) ; 
    if ( ! ofs.is_open() ) {
      std::cerr << "Cannot open file." << std::endl ;
      return -1 ;
    }
  }
  // polyNum, total consNum, redundancy, variNum, zeros
   
  double totalTime = 0 ;
  int idx = 0 ;
  std::string polyName = "P_" + std::to_string( ioInter.get_cons_num() )
      + "_" + std::to_string( ioInter.get_redundancy() )
      + "_" + std::to_string( ioInter.get_vari_num() )
      + "_" + std::to_string( ioInter.get_zero_num() ) ;
  
  std::cout << polyName << std::endl ;
  for (auto& poly : polyVec) {
    if (argc == 3) {
      ofs << ioInter.get_cons_num() << ","
        << poly.get_redundant_num() << "," << ioInter.get_vari_num() << ","
        << ioInter.get_zero_num() << "," ; 
    }
    std::cout << "poly: " << idx << std::endl ;
    auto start = std::chrono::steady_clock::now() ;
    poly.MinimizeSimple() ;
    auto end = std::chrono::steady_clock::now() ;
    std::chrono::duration<double> diff = (end - start) * 1000 ;
    std::cout << "Total running time: " << diff.count() << " ms" << std::endl ;
    if (argc == 3) {
      ofs <<  idx << ", " << diff.count() << std::endl ;
    }
    ++ idx ;
    totalTime += diff.count() ;
    poly.PrintActiveIdx() ;
  }
  double aveTime = totalTime / polyVec.size() ;
  std::cout << "Average running time: " << aveTime << " ms" << std::endl ;
  if (argc == 3) {
    ofs << aveTime << std::endl ; 
  }

  return 0 ;
} 
