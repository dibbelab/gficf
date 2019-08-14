#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppEigen.h>

// Compute jaccard coefficient between nearest-neighbor sets in parallell
//
// Author: Gennaro Gambardella, Date: 12/08/2019

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
struct JCoefficient : public RcppParallel::Worker {
  
  // input matrix to read from
  const Rcpp::NumericMatrix mat;
  
  // output matrix to write to
  Eigen::MatrixXd& rmat;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  JCoefficient(const Rcpp::NumericMatrix& mat, Eigen::MatrixXd& rmat)
    : mat(mat), rmat(rmat) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < mat.ncol(); j++) {
        
        int k = mat(i,j)-1;
        
        Rcpp::NumericVector nodei = mat(i,Rcpp::_);
        Rcpp::NumericVector nodej = mat(k,Rcpp::_);
        int u = intersect(nodei, nodej).size();  // count intersection number

        if(u>0){
          rmat((i*mat.ncol())+j, 0) = i+1;
          rmat((i*mat.ncol())+j, 1) = k+1;
          rmat((i*mat.ncol())+j, 2) = u/(2.0*mat.ncol() - u);
        }
      }
    }
  }
};

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_parallel_jaccard_coef(Rcpp::NumericMatrix mat, bool printOutput) {
  
  if (printOutput)
  {
    Rprintf("Running Parallell Jaccard Coefficient Estimation...\n");
  }
  // allocate the matrix we will return
  Eigen::MatrixXd rmat = Eigen::MatrixXd::Zero(mat.nrow()*mat.ncol(),3);

  // create the worker
  JCoefficient JCoefficient(mat, rmat);
  
  // call it with parallelFor
  parallelFor(0, mat.nrow(), JCoefficient);
  
  if (printOutput)
  {
    Rprintf("Done!!\n");
  }
  return rmat;
}