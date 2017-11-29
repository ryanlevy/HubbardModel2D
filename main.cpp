#define DEBUG 1

#include <iostream>
#include <iomanip>
#include <string>
//#include <cmath>
#include "Hubbard2D.hpp"

#ifdef USE_SPECTRA
#include <MatOp/SparseGenMatProd.h>
#include <SymEigsSolver.h>
#endif

#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <ietl/interface/eigen3.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>

///check if a file exists
bool fileExists(const std::string& fname){
  std::ifstream f(fname.c_str());
  bool isGood = f.good();
  f.close();
  return isGood;
}

int main(int argc, const char * argv[]) {
  int ne; //# of electrons
  double U;

  if(argc == 1){
    std::cerr<<"Using all defaults"<<std::endl;
    ne=2; //# of electrons
    U=4;
  }
  else if(argc >5 || (argc-1)%2){
    std::cout<<"Wrong input format!"<<std::endl;
    std::cout<<"-ne [#num electrons] -U [# for U/t]"<<std::endl;
    return 1;
  }
  else{
    ne=2;
    U=4;
    for(int i=1;i<argc;i+=2){
     if(std::string(argv[i]) == "-ne")
       ne = std::atoi(argv[i+1]);
     if(std::string(argv[i]) == "-U")
       U = std::atof(argv[i+1]);
    }
  }
  std::cout << "2D "<<BASISSIZE/2<<" site Hubbard model with "
  << ne <<" up/down electrons and U/t="<<U<<std::endl;


  HubbardModel2D H(ne,ne,1,U);
  H.makeBasis();
  std::cout <<"Made Basis: " << H.getBasis()->size()<< std::endl;

  H.buildHubbard2D();
  std::cout <<"Built Matrix" << std::endl;
  SpMat Hmat = *(H.getH());
  std::cout << "# of non-zero elements:"<<(H.getH())->nonZeros() << std::endl;
  H.getH()->makeCompressed();

  std::ofstream outfile;
  if(fileExists("energies.dat")){
    outfile.open("energies.dat",std::ofstream::out |std::ofstream::app);
  }
  else{
     outfile.open("energies.dat",std::ofstream::out |std::ofstream::app);
     outfile<<"#ne U E0"<<std::endl;
  }
  outfile << std::setprecision(9);

  std::ofstream outfile_v;
  if(fileExists("ketvals.dat")){
    outfile_v.open("ketvals.dat",std::ofstream::out |std::ofstream::app);
  }
  else{
    outfile_v.open("ketvals.dat",std::ofstream::out |std::ofstream::app);
    outfile_v <<"#ne U dbl up|down antiferro"<<std::endl;
  }
  outfile_v << std::setprecision(9);

  //construct test states
  size_t spinOffset = BASISSIZE/2;

  lattice_t psi1;
  int i;
  for(i=0;i<ne;i++)
    psi1[i].flip();

  for(i=0;i<ne;i++)
    psi1[(i+spinOffset)%(BASISSIZE)].flip();

  lattice_t psi2;
  i=0;

  for(i=0;i<ne;i++)
    psi2[i].flip();

  for(i=0;i<ne;i++)
    psi2[(i+ne+spinOffset)%(BASISSIZE)].flip();
  
  lattice_t psi3;
  i=0;
  for(i=0;i<ne;i++)
    psi3[2*i].flip();
  for(i=0;i<ne;i++)
    psi3[(2*i+1+spinOffset)%(BASISSIZE)].flip();
  //std::cout<<psi1<<std::endl;
  //std::cout<<psi2<<std::endl;
  //std::cout<<psi3<<std::endl;

  std::string double_occ="|";
  for(i=0;i<ne;i++)
    double_occ+="↑↓|";
  std::string not_occ="|";
  for(i=0;i<ne;i++)
    not_occ+="↑|";
  for(i=0;i<ne;i++)
    not_occ+="↓|";
  std::string antiferro="|";
  for(i=0;i<ne;i++)
    antiferro+="↑|↓|";

  //**************************************
  //        generate eigenvalues
  //**************************************

  //Diag with specrta
#ifdef USE_SPECTRA
  // Construct matrix operation object using the wrapper class SparseGenMatProd
  using wrapper_t =Spectra::SparseGenMatProd<double,Eigen::RowMajor,long>;
  wrapper_t op(Hmat);
  
  // Construct eigen solver object, requesting the smallest eigenvalue
  Spectra::SymEigsSolver< double, Spectra::SMALLEST_ALGE , wrapper_t > eigs(&op, 1, 6);
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
#ifdef DEBUG
  if(eigs.info() != Spectra::SUCCESSFUL)
      std::cout <<"Warning something failed " <<nconv << std::endl;
#endif

  std::cout<< "E0="<<eigs.eigenvalues()[0]<<std::endl;

  std::cout<<double_occ<<" :"<<eigs.eigenvectors()(H.getState(psi1),0)<<std::endl; //psi1
  std::cout<<not_occ<<" :"<<eigs.eigenvectors()(H.getState(psi2),0)<<std::endl; //psi2
#else
  typedef boost::lagged_fibonacci607 Gen;
  Gen mygen;
  mygen.seed(0);
  
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor,long> Matrix;
  typedef Eigen::VectorXd Vector;
  typedef ietl::vectorspace<Vector> Vecspace;
   
   
   
  // Creation of an iteration object:
  int max_iter = 10*H.getBasis()->size();
  double rel_tol = 500*std::numeric_limits<double>::epsilon();
  double abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3);
  int n_lowest_eigenval = 1;
  std::vector<double> eigen;
  std::vector<double> err;
  std::vector<int> multiplicity;

  std::vector<double> groundStates;

  Vecspace vec(Hmat.cols());
  ietl::lanczos<Matrix,Vecspace> lanczos(Hmat,vec);
  ietl::lanczos_iteration_nlowest<double>
  iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
  std::cout << "Running lanczos" << std::endl;
  try{
    lanczos.calculate_eigenvalues(iter,mygen);
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
    std::cout<<"number of iterations: "<<iter.iterations()<<"\n";
    groundStates.push_back(eigen[0]);
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << "\n";
  }
  std::cout << "#        eigenvalue            error         multiplicity\n";
  for (int i=0;i<10;++i)
    std::cout << i+1 << "\t" << eigen[i] << "\t" << err[i] << "\t"
    << multiplicity[i] << "\n";

  outfile<<ne<<" " <<U<< " "<< groundStates[0]<<std::endl;
  //**************************************
  //        generate eigenvectors
  //**************************************
  // call of eigenvectors function follows:
  std::cout << "\nHead of Eigenvector for the lowest eigenvalue:\n\n";
  std::vector<double>::iterator start = eigen.begin();
  std::vector<double>::iterator end = eigen.begin()+1;
  std::vector<Vector> eigenvectors; // for storing the eigen vectors.
  ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).

  try {
    lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen,max_iter);
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << "\n";
  }
  std::cout << eigenvectors[0].head(10)<<std::endl;
  std::cout << "Information about the eigenvector computations:\n\n";
  for(int i = 0; i < info.size(); i++) {
    std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
    << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
    << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
    << info.residual(i) << " error_info(" << i+1 << "): "
    << info.error_info(i) << "\n\n";
  }

  std::cout<<double_occ<<" :"<<eigenvectors[0](H.getState(psi1))<<std::endl; //psi1
  std::cout<<not_occ<<" :"<<eigenvectors[0](H.getState(psi2))<<std::endl; //psi2
  std::cout<<antiferro<<" :"<<eigenvectors[0](H.getState(psi3))<<std::endl; //psi3

  outfile_v << ne << " "<<U<<" " << eigenvectors[0](H.getState(psi1)) << " " 
            <<eigenvectors[0](H.getState(psi2))<<" " <<eigenvectors[0](H.getState(psi3))<<std::endl;
#endif

  return 0;
}
