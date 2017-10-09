//
//  Hubbard2D.hpp
//  HubbardModel
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

#define BASISSIZE 32 //8 //2*4
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
    // TODO:allow for t' values
    Eigen::MatrixXd HH;

    std::vector<std::vector<int> > neighbors; 
    
public:
    HubbardModel2D(int Lx_in,int Ly_in,int nUp_in,int nDown_in,double t_in, double U_in):
    Lx(Lx_in),Ly(Ly_in),nUp(nUp_in),nDown(nDown_in),t(t_in),U(U_in){
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

    unsigned long getState(const lattice_t &psi){return indexMap[psi];}
    
};

#endif /* Hubbard2D_hpp */
