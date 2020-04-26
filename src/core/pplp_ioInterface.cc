/*******************************************************************************
 * Copyright (c) 2016 Dec. Verimag. All rights reserved.
 * @author Hang YU
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <regex>
#include <algorithm>
#include "pplp_ioInterface.h"

using namespace PPLP ;

IoInterface::IoInterface() 
  : _cons_num(0), _vari_num(0), _zero_num(0)
{}

/*******************************************************************************
 * Read polyhedron from a file. The polyhedron is stored as matrix in that file.
*******************************************************************************/
std::vector<Polyhedron> IoInterface::LoadPolyhedra(const char* filepath) {
  std::vector<Polyhedron> polyVec ;
#ifndef NO_STRING
  std::ifstream ifs(filepath);
  if ( ! ifs.is_open() ) {
    std::cerr << "LoadPolyhedra: Cannot open file " 
        << std::string(filepath) << std::endl ;
    std::terminate() ;
  }
  std::string line ;
  std::string infoName ;
  int infoVal ;
  int pos ;
  std::regex pattern("[[:space:]]") ;
  while( std::getline(ifs, line) ) {
    if (line.length() == 0) continue ;
    line = std::regex_replace(line, pattern, "");
    if (line == "begin" || line == "Begin" || line == "BEGIN") break ;
    pos = line.find("=") ;
    infoName = line.substr(0, pos) ;
    infoVal = std::stoi(line.substr(pos+1));
    // tolower is a global function
    std::transform(infoName.begin(), infoName.end(), infoName.begin(), ::tolower);
    if (infoName == "constraints") {
      _cons_num = infoVal ;
    }
    if (infoName == "variables") {
      _vari_num = infoVal ;
    }
    if (infoName == "zeros") {
      _zero_num = infoVal ;
    }
    if (infoName == "redundancy") {
      _redundancy = infoVal ;
    }
  }

  if (_cons_num == 0 || _vari_num == 0) {
    std::cerr << "Need information of the number of Polyhedra, constraints and variables" << std::endl ;
    std::terminate() ;
  }

  while ( ! ifs.eof() ) {

    Matrix coeff(_cons_num, _vari_num+1) ;
    std::vector<int> ineqIdx ;
    std::vector<int> eqIdx ;

    // read the name of polyhedron
    while( ! ifs.eof() && std::getline(ifs, line) ) {
      line = std::regex_replace(line, pattern, "") ;
      if (line[0] >= 'a' && line[0] <= 'z') break ;
      else if (line[0] >= 'A' && line[0] <= 'Z') break ;
      else if (line[0] >= '0' && line[0] <= '9') {
        std::cerr << "The name of Polyhedron is missing or illegal." << std::endl ;
        std::terminate() ;
      }
      else if (line.length() != 0) {
        std::cerr << "The name of Polyhedron is illegal." << std::endl ;
        std::terminate() ;
      }
    }
    if ( ifs.eof() ) break ;
    // the polyhedra names are in the form: P_C_R_V_Z_G_I, 
    // where C: total number of constraints,
    // R: number of redundant constraints, V: number of variables
    // G: number of generators, I: the identity of polydedron
    //Polyhedron newPoly(_cons_num, _vari_num) ;
    int genNum, id ;
    std::regex rgx ("_") ; 
    // default constructor: end of sequence
    std::regex_token_iterator<std::string::iterator> end ; 
    std::regex_token_iterator<std::string::iterator> 
      it (line.begin(), line.end(), rgx, -1) ;
    // ignore the information except the number of generator and id
    // TODO: maybe do not need the information of ioInterface class
    for (int p = 0; p < 5; ++ p) {
      ++ it ;
    }
    //newPoly.set_generator_num( std::stoi(*it) ) ;
    genNum = std::stoi(*it) ;
    ++ it ;
    //newPoly.set_id( std::stoi(*it) ) ;
    id = std::stoi(*it) ;
    // don't set redundant number, but use minimize() to calculate
    //newPoly.set_redundant_num(0) ;
    //newPoly.set_zero_num(_zero_num) ;
    // read the constraints matrix
    double newVal ;
    for (int i = 0; i < _cons_num; i++) {
      char op = ' ' ;
      for (int j = 0; j < _vari_num; j++) {
        ifs >> newVal ;
        //newPoly.SetCoef(i, j, newVal) ;
        coeff(i, j) = newVal ;
      }
      // read the operators begin
      while (op == ' ') {
        ifs >> op ;
      }
      if (op == '=') {
        //newPoly.SetOperator(i, ConstraintOperator::equal) ; 
        eqIdx.push_back(i) ;
      }
      else if (op == '<') {
        ifs >> op ;
        if (op == '=') {
          //newPoly.SetOperator(i, ConstraintOperator::lesseq) ;
          ineqIdx.push_back(i) ;
        }
        /*
        // there is no < currently
        else {
          newPoly.SetOperator(i, ConstraintOperator::less) ;
        }
        */
      } 
      // read the operators end
      ifs >> newVal ;
      //newPoly.SetConstant(i, newVal) ;
      coeff(i, _vari_num) = newVal ;
    }
    // divide the coefficients by their gcd
    double curr ;
    for (int i = 0 ; i < _cons_num; ++ i) {
      int gcd = Tool::GetGcd( coeff.row(i).head(_vari_num), coeff(i, _vari_num) ) ;
      if (gcd == 1) continue ;
      for (int j = 0; j < _vari_num+1; ++ j) {
        curr = coeff(i, j) / gcd ;
        coeff(i, j) = curr ;
      }
    }
    
#ifdef DEBUGINFO_RAYTRACING
  std::cout << "Load new polyhedron: " << std::endl ;
  newPoly.Print() ;
#endif

    int eqConsNum = eqIdx.size() ;
    int ineqConsNum = ineqIdx.size() ;
    Matrix eqMatrix(eqConsNum, _vari_num+1) ;
    Matrix ineqMatrix(ineqConsNum, _vari_num+1) ;
    for (int i = 0; i < eqConsNum; ++ i) {
      for (int j = 0; j < _vari_num+1; ++ j) {
        eqMatrix(i, j) = coeff(eqIdx[i], j) ;
      }
    }
    for (int i = 0; i < ineqConsNum; ++ i) {
      for (int j = 0; j < _vari_num+1; ++ j) {
        ineqMatrix(i, j) = coeff(ineqIdx[i], j) ;
      }
    }


    int nonDupIneqNum, nonDupEqNum ;
    std::vector<int> nonDupIneqIdx, nonDupEqIdx ;
    if (ineqConsNum != 0) {
      nonDupIneqIdx = Tool::GetNonDupIdx(ineqMatrix) ;
      nonDupIneqNum = nonDupIneqIdx.size() ;
    }
    else {
      nonDupIneqNum = 0 ;
    }
    if (eqConsNum != 0) {
      nonDupEqIdx = Tool::GetNonDupIdx(eqMatrix) ; 
      nonDupEqNum = nonDupEqIdx.size() ;
    }
    else {
      nonDupEqNum = 0 ;
    }

    Polyhedron newPoly(nonDupIneqNum, _vari_num, nonDupEqNum) ;

    for (int i = 0; i < nonDupIneqNum; ++ i) {
      for (int j = 0; j < _vari_num; ++ j) {
        curr = ineqMatrix(nonDupIneqIdx[i], j) ;
        newPoly.SetCoef(i, j, curr) ;
      }
      curr = ineqMatrix(nonDupIneqIdx[i], _vari_num) ;
      newPoly.SetConstant(i, curr) ;
    }
    for (int i = 0; i < nonDupEqNum; ++ i) {
      for (int j = 0; j < _vari_num; ++ j) {
        curr = eqMatrix(nonDupEqIdx[i], j) ;
        newPoly.SetEqCoef(i, j, curr) ;
      }
      curr = eqMatrix(nonDupEqIdx[i], _vari_num) ;
      newPoly.SetEqConstant(i, curr) ;
    }

    newPoly.set_generator_num(genNum) ;
    newPoly.set_zero_num(_zero_num) ;
    newPoly.set_redundant_num(0) ;
    newPoly.set_id(id) ;

#ifdef DEBUGINFO_RAYTRACING
  std::cout << "After removing duplicated constraints: " << std::endl ;
  newPoly.Print() ;
#endif

    polyVec.push_back( std::move(newPoly) ) ;
  }
#endif
  return polyVec ;
}

int IoInterface::get_cons_num() {
  return _cons_num ;
}

int IoInterface::get_vari_num() {
  return _vari_num ;
}

int IoInterface::get_zero_num() {
  return _zero_num ;
}

int IoInterface::get_redundancy() {
  return _redundancy ;
}

/*
void IoInterface::WritePoly(const std::string& filepath, const Polyhedron& poly) {
#ifndef NO_STRING
  int id = 0 ;
  std::ifstream ifs(filepath) ;
  bool fileExists = true ;
  // the file does not exist
  if ( ifs.fail() ) {
    fileExists = false ;
  }
  else {
    // read the file to get id
    std::string line ;
    std::string infoName ;
    int inConsNum, pos, infoVal ;
    std::regex pattern("[[:space:]]") ;
    while( std::getline(ifs, line) ) {
      if (line.length() == 0) continue ;
      line = std::regex_replace(line, pattern, "");
      if (line == "begin" || line == "Begin" || line == "BEGIN") break ;
      pos = line.find("=") ;
      infoName = line.substr(0, pos) ;
      infoVal = std::stoi(line.substr(pos+1));
      // tolower is a global function
      std::transform(infoName.begin(), infoName.end(), infoName.begin(), ::tolower);
      if (infoName == "constraints") {
        inConsNum = infoVal ;
      }
    }
    while ( ! ifs.eof() ) {
      // read the name of polyhedron
      while( ! ifs.eof() && std::getline(ifs, line) ) {
        line = std::regex_replace(line, pattern, "") ;
        if (line[0] >= 'a' && line[0] <= 'z') {
          ++ id ;
          break ;
        }
        else if (line[0] >= 'A' && line[0] <= 'Z') {
          ++ id ;
          break ;
        }
        else if (line[0] >= '0' && line[0] <= '9') {
          std::cerr << "The name of Polyhedron is missing or illegal." << std::endl ;
          std::terminate() ;
        }
        else if (line.length() != 0) {
          std::cerr << "The name of Polyhedron is illegal." << std::endl ;
          std::terminate() ;
        }
      }
      if ( ifs.eof() ) break ;
      for (int i = 0; i < inConsNum; ++ i) {
        std::getline(ifs, line) ;
      } 
    }
    ifs.close() ;
  }
  std::ofstream ofs(filepath, std::ofstream::out | std::ofstream::app) ; 
  if ( ! ofs.is_open() ) {
    std::cerr << "WritePoly: Cannot open file " 
        << std::string(filepath) << std::endl ;
    std::terminate() ;
  }
  int consNum = poly.get_constraint_num() ;
  int reConsNum = poly.get_redundant_num() ;
  int variNum = poly.get_variable_num() ;
  int zeroNum = poly.get_zero_num() ;
  if (fileExists == false) {
    std::string line ;
    line = "constraints = " + std::to_string(consNum) ;
    ofs << line << std::endl ;
    line = "redundancy = " + std::to_string(reConsNum) ;
    ofs << line << std::endl ;
    line = "variables = " + std::to_string(variNum) ;
    ofs << line << std::endl ; 
    line = "zeros = " + std::to_string(zeroNum) ;
    ofs << line << std::endl ;  
    ofs << std:: endl ;
    ofs << "Begin" << std::endl << std:: endl ;
  }
// write polyhedra name
  std::string name = "P_" + std::to_string(consNum) 
    + "_" + std::to_string(reConsNum) 
    + "_" + std::to_string(variNum) 
    + "_" + std::to_string(zeroNum) 
    + "_" + std::to_string( poly.get_generator_num() ) 
    + "_" + std::to_string(id) ;
  ofs << name << std::endl ;
  for (int i = 0; i < consNum; ++ i) {
    for (int j = 0; j < variNum; ++ j) {
      ofs << poly.GetCoef(i, j) << " " ;
    }
    // TODO consider different operators later
    ofs << "<=" << " " << poly.GetConstant(i) << std::endl ;
  }
  ofs << std::endl ;
  ofs.close() ;
#endif
}
*/

