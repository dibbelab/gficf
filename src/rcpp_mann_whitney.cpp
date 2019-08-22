#include "mann_whitney.h"
#include <Rcpp.h>

// Compute two-sided Mann–Whitney U test between two umpaired samples. 
// In R this c++ function corresond to wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = F)
//
// Author: Gennaro Gambardella, Date: 20/08/2019



// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_WMU_test(Rcpp::NumericMatrix M,Rcpp::NumericVector idx1,Rcpp::NumericVector idx2) {
  Rcpp::NumericMatrix res(M.nrow(),2);
  
  for(int i=0;i<M.nrow();i++)
  {
    Rcpp::NumericVector g1 = subset(M.row(i),idx1);
    Rcpp::NumericVector g2 = subset(M.row(i),idx2);
    res(i,0) = MWUtest(g1,g2);
    res(i,1) = std::log2( (double) (avg(g1+1)/avg(g2+1)) );
  }
  
  return res;
}
