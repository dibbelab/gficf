#include "mann_whitney.h"
#include <cmath>
#include <Rcpp.h>
#include <RcppParallel.h>

// Compute two-sided Mann–Whitney U test with coninuity correction between two umpaired samples.
// In R this c++ function corresond to wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = F)
//
// Author: Gennaro Gambardella, Date: 20/08/2019

// [[Rcpp::depends(RcppParallel)]]
struct WMU_test : public RcppParallel::Worker {
  
  // input matrix to read from
  const RcppParallel::RMatrix<double> matX;
  const RcppParallel::RMatrix<double> matY;

  // output matrix to write to
  RcppParallel::RMatrix<double> rmat;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  WMU_test(const Rcpp::NumericMatrix& matX, const Rcpp::NumericMatrix& matY, Rcpp::NumericMatrix& rmat)
    : matX(matX), matY(matY), rmat(rmat) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; k++) {
        
        RcppParallel::RMatrix<double>::Row row1(matX.row(k));
        RcppParallel::RMatrix<double>::Row row2(matY.row(k));
      
        std::vector<double> v1(row1.length());
        std::copy(row1.begin(), row1.end(), v1.begin());
        
        std::vector<double> v2(row2.length());
        std::copy(row2.begin(), row2.end(), v2.begin());
        
        double pval = 1;
        std::vector<double> values(v1.size());
        std::vector<double> valuesOrd(v1.size()+v2.size());
        std::vector<double> valuesRnk(v1.size()+v2.size());
        
        // concat vectors
        std::copy(v1.begin(), v1.end(), values.begin());
        values.insert( values.end(),v2.begin(),v2.end());
        
        // sort
        std::vector<size_t> idx = sort_indexes(values);
        for(size_t i=0;i<values.size();i++)
        {
          //std::cout << values[idx[i]] << std::endl;
          valuesOrd[i] = values[idx[i]];
        }
        // get ranks
        valuesRnk = getRanks(valuesOrd);
        
        std::vector<double> nties = getCounts(valuesOrd);
        
        if(nties.size()>1) {
          
          double long U1 = (v1.size() * (v1.size()+1)) * -0.5; 
          double long U2 = (v2.size() * (v2.size()+1)) * -0.5;
          std::vector<double> u1_rnk(v1.size(),0);
          std::vector<double> u2_rnk(v2.size(),0);
          std::vector<double> v1_ord(v1.size(),0);
          std::vector<double> v2_ord(v2.size(),0);
          int u1_i=0;
          int u2_i = 0;
          
          for(size_t i=0;i<valuesOrd.size();i++)
          {
            if(idx[i]<v1.size())
            {
              U1 +=valuesRnk[i];
              u1_rnk[u1_i] = valuesRnk[i];
              v1_ord[u1_i] = valuesOrd[i];
              u1_i++;
            } else {
              U2 +=valuesRnk[i];
              u2_rnk[u2_i] = valuesRnk[i];
              v2_ord[u2_i] = valuesOrd[i];
              u2_i++;
            }
          }
          
          double mu = v1.size()*v2.size()/2, sig=0, z = 0, correction = 0;
          sig = getSigma(nties,v1.size(),v2.size());
          
          z = U1<U2 ? U1-mu : U2-mu;
          z = z<0 ? z+0.5 : z-0.5; // Continuity Correction
          z = z/sig;
          
          pval = getPvalue(z);
        }
        
        rmat(k,0) = pval;
        transform(v1.begin(), v1.end(), v1.begin(),bind2nd(std::plus<double>(), 1.0));
        transform(v2.begin(), v2.end(), v2.begin(),bind2nd(std::plus<double>(), 1.0));
        rmat(k,1) = log2( (double) (avg(v1)/avg(v2)) );
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_parallel_WMU_test(Rcpp::NumericMatrix matX,Rcpp::NumericMatrix matY, bool printOutput) {
  
  if (printOutput)
  {
    Rprintf("Running Parallell WM-U test...\n");
  }
  
  // allocate the matrix we will return
  Rcpp::NumericMatrix rmat(matX.nrow(),2);
  
  // create the worker
  WMU_test WMU_test(matX, matY, rmat);
  
  // call it with parallelFor
  parallelFor(0, matX.nrow(), WMU_test);
  
  if (printOutput)
  {
    Rprintf("Done!!\n");
  }
  return rmat;
}
