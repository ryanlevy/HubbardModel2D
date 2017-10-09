//
//  Hubbard2D.cpp
//
//  Created by ryle on 1/15/17.
//  Copyright © 2017 R L. All rights reserved.
//

#include <cassert>
#include "Hubbard2D.hpp"
#include <iostream>

unsigned long nextBitPerm(unsigned long v){
    //https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
    //unsigned int z = (v | (v - 1)) + 1; // t gets v's least significant 0 bits set to 1
    unsigned long z = v | (v - 1);
    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    //return v==0 ? 0 : z | ((((z & -z) / (v & -v)) >> 1) - 1);
    return (z + 1) | (((~z & -~z) - 1) >> (__builtin_ctzl(v) + 1));
}
void HubbardModel2D::init(){
  std::cerr<<"INIT: HAMILTONIAN HUBBARD "<<std::endl;
  neighbors.resize(BASISSIZE/2);
  std::ifstream infile;
  infile.open("bond.dat");
  assert(infile);
  while (!infile.eof()){
    int i; int j;
    infile>>i;
    infile>>j;
    if (!infile.eof()){
      neighbors[i].push_back(j);
      //neighbors[j].push_back(i);
    }
  }
  infile.close();    

  std::string vLoc="vec.dat";
  //std::string vLoc="";

  int L =BASISSIZE/2; //trust defs.h will be setup properly
  HH.resize(L,L);
  HH.setZero();
  
  //build the hopping matrix
  for(int i=0;i<neighbors.size();i++){
    auto nlist = neighbors[i];
    for(auto n = nlist.begin();n!=nlist.end();n++){
      HH(i,*n)=-t;
    }
  } 
  //T is hermitian
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      assert( std::abs(HH(i,j)-HH(j,i)) < 1e-6);
    }
  } 
#ifdef DEBUG
  std::cout << "T matrix is Hermitian!" << std::endl;
#endif

  std::cerr<<"DONE "<<std::endl;

}
void HubbardModel2D::makeBasis(){
    assert(nUp+nDown <=BASISSIZE); //sanity check
    
    //uneffecient but good enough for now
    
    using config_t = unsigned long;
    
    std::vector<config_t > upConfigs;
    std::vector<config_t > downConfigs;
    
    //start with ...000[1]*nUp, smallest valid config
    config_t  startUp = (1<<nUp) -1;
    config_t  vup=0;
    config_t  nextup = startUp;
    
    while(nextup <(1<<(BASISSIZE/2))){
        vup = nextup;
        nextup = nextBitPerm(vup);
        upConfigs.push_back(vup);
    }
    
    config_t  startDown = (1<<nDown) -1;
    config_t  vDown=0;
    config_t  nextDown = startDown;
    
    while(nextDown <(1<<(BASISSIZE/2))){
        vDown = nextDown;
        nextDown = nextBitPerm(vDown);
        downConfigs.push_back(vDown);
    }

    //candidates assembled, build basis
    using Iter = std::vector<config_t >::const_iterator;

    for(Iter up=upConfigs.begin(); up!=upConfigs.end(); ++up){
        for(Iter down=downConfigs.begin(); down!=downConfigs.end(); ++down){
            //[up][down]
            lattice_t config((*up<<(BASISSIZE/2))+*down);
            basis.push_back(config);
        }
    }//end basis loop
    //std::cout << basis.size() << std::endl;
    //std::cout << basis[0]<< std::endl;
}

/// count the electrons between start and end; offset =1 if down 0 if up
int countElectrons(const lattice_t &psi,size_t start,size_t end,size_t offset){
  int count=0;
  if(start>end)
    std::swap(start,end);

  for(size_t i=start;i<=end;i+=1){
    if(psi[i+offset])
      count++;
  }
  //we counted ourselves, need to subtract 1
#ifdef DEBUG
  assert(((psi[end+offset]==1) ^ (psi[start+offset]==1)) || start==end);
  assert(count>0);
#endif
  return count-1;
}

void HubbardModel2D::buildHubbard2D(){
    
    double epsilon =0; //1e-30;

    typedef Eigen::Triplet<double, unsigned long> T;
    std::vector<T> tripletList;
    
    size_t nH = basis.size();
    assert(nH!=0);

    H.resize(nH,nH);
    tripletList.reserve((unsigned int) (nH*0.05));
    
    
    using Iter = basis_t::const_iterator;
    if(indexMap.size() == 0){
        indexMap.clear();
        
        size_t i=0;
        for(Iter it=basis.begin(); it!=basis.end();it++,i++)
            indexMap[*it]=i;
    }
    
    int spinOffset=BASISSIZE/2;
    // Hit each |\psi> with H
#pragma omp declare reduction (merge : std::vector<T> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for reduction(merge: tripletList)
    for(size_t basisIdx=0;basisIdx<nH;basisIdx++){
    //for(Iter it=basis.begin();it!=basis.end();it++){
        lattice_t psi=basis[basisIdx];//*it;
        unsigned long ci = indexMap[psi];
        double diagWeight = 0.0;
        //first do the hopping matrix, which need to know about
        //neighbors though T
        for(int loc=0;loc<BASISSIZE/2;loc++){
          // each site is (up,down)
          //spin up then spin down
          for(size_t offset : {0,spinOffset}){
            if(psi[loc+offset]==1){
              //we have a spin here, now we need to go through all places it can hop
              for(int Li=0;Li<BASISSIZE/2;Li++){
                double hopVal = HH(loc,Li); //probably need to be more careful if complex
                bool occupied = psi[Li+offset]==1;
                if(abs(hopVal)>epsilon && (!occupied || Li==loc)){
                  if(Li==loc){
                    //assert(false); //unsure that this should ever come up in a hopping matrix
                    diagWeight+=hopVal;
                  }
                  else{
                    lattice_t phi=psi;
                    phi[loc+offset].flip();
                    phi[Li+offset].flip();
                    //sign due to normal ordering
                    int sign = (countElectrons(psi,loc,Li,offset)%2 ==0) ? 1 : -1;
                    tripletList.push_back(T(ci, indexMap[phi], sign*hopVal));
                  }
                }
              }//end lattice hop
            }//end occupied
          }//end offset
       
          //U Term
          if(psi[loc]==1 && psi[loc+spinOffset]==1)
            diagWeight+=U;
        }//linear site loop

        tripletList.push_back(T(ci,ci,diagWeight));

    }//end H|\psi>
    H.setFromTriplets(tripletList.begin(), tripletList.end());
}

