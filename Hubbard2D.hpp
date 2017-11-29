//
//  Hubbard2D.hpp
//
//  Created by ryle on 1/15/17.
//  Copyright © 2017 R L. All rights reserved.
//

#ifndef Hubbard2D_hpp
#define Hubbard2D_hpp

#include <stdio.h>
#include <bitset>
#include <complex>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <Eigen/SparseCore>
#include <Eigen/Core>
//**************************************c
#define BASISSIZE 32 //(4*4)*2
//**************************************c
typedef std::bitset<BASISSIZE> tbitset;
typedef tbitset lattice_t;
typedef std::vector<lattice_t> basis_t;
typedef std::complex<double> Complex;

typedef Eigen::SparseMatrix<double,Eigen::RowMajor,long> SpMat;

///basis mapping
typedef std::unordered_map<lattice_t,unsigned long> dict_t;
typedef Eigen::MatrixXd rotation_t;

EIGEN_STATIC_ASSERT(BASISSIZE%2==0, BASISSIZE_NEEDS_TO_BE_EVEN);

class HubbardModel2D{
protected:
    int Lx,Ly,nUp,nDown,numParam;
    double t,U;
    basis_t basis;
    
    //Hamiltonian
    SpMat H;
    ///make a map of basis bitset to index
    dict_t indexMap;

    /// Hoping/Neighbor Matrix; should include -t,-t' etc
    // TODO:allow for t' values that aren't 1
    Eigen::MatrixXd HH;

    //stores neighbors, currently unused beyond HH generation
    std::vector<std::vector<int> > neighbors; 
    
public:
    HubbardModel2D(int nUp_in,int nDown_in,double t_in, double U_in):
    nUp(nUp_in),nDown(nDown_in),t(t_in),U(U_in){
      init();};
    
    //setup matrices
    void init();
    ///construct a basis given Lx,Ly, nUp, nDown
    void makeBasis();
    ///build matrix given init bonds
    void buildHubbard2D();
    
    //getter functions
    
    ///getter for the internal basis
    const basis_t *getBasis() {return &basis;}
    ///getter for last made Hamiltonian
    SpMat *getH() {return &H;}

    ///given a basis element return index
    unsigned long getState(const lattice_t &psi){return indexMap[psi];}
    
};

#endif /* Hubbard2D_hpp */
