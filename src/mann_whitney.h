#include <Rcpp.h>

// Compute two-sided Mann–Whitney U test between two umpaired samples. 
// In R this c++ function corresond to wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = F)
//
// Author: Gennaro Gambardella, Date: 20/08/2019


template <typename T>
std::vector<size_t> sort_indexes(std::vector<T> const &v);

std::vector<double> getRanks(std::vector<double> const &absoluteValues);

std::vector<double> getUniq(std::vector<double> u);

std::vector<double> getCounts(std::vector<double> const &valuesOrd);

double getSigma(std::vector<double>& u_counts,double n1, double n2);

double getPvalue(double& z);

Rcpp::NumericVector subset(Rcpp::NumericVector const &v, Rcpp::NumericVector idx);

double avg(Rcpp::NumericVector const &v);

double avg(std::vector<double> const &v);

double MWUtest(Rcpp::NumericVector const &v1, Rcpp::NumericVector const &v2);
