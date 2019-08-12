#include <Rcpp.h>
#include <RcppParallel.h>

// Compute jaccard coefficient between nearest-neighbor sets in parallell
//
// Author: Gennaro Gambardella, Date: 12/08/2019


// [[Rcpp::depends(RcppParallel)]]
struct JCoefficient : public RcppParallel::Worker {
  
  // input matrix to read from
  const RcppParallel::RMatrix<double> mat;
  
  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  JCoefficient(const Rcpp::NumericMatrix& mat, Rcpp::NumericMatrix& rmat)
    : mat(mat), rmat(rmat) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < mat.ncol(); j++) {
        
        int k = mat(i,j)-1;
        std::vector<double> v(mat.ncol());
        std::vector<double>::iterator it;
        
        RcppParallel::RMatrix<double>::Row v1 = mat.row(i);
        RcppParallel::RMatrix<double>::Row v2 = mat.row(k);
        
        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());
        
        std::vector<int> v_intersection;
        
        std::set_intersection(v1.begin(), v1.end(),
                              v2.begin(), v2.end(),
                              std::back_inserter(v_intersection));
        int u = v_intersection.size();

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
Rcpp::NumericMatrix rcpp_parallel_jaccard_coef(Rcpp::NumericMatrix mat, bool printOutput) {
  
  if (printOutput)
  {
    Rprintf("Running Parallell Jaccard Coefficient Estimation...\n");
  }
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(mat.nrow()*mat.ncol(),3);
  
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