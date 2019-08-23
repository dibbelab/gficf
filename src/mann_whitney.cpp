#include "mann_whitney.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_cdf.h>

// Compute two-sided Mann–Whitney U test with coninuity correction between two umpaired samples. 
// In R this c++ function corresond to wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = T)
//
// Author: Gennaro Gambardella, Date: 20/08/2019


template <typename T>
std::vector<size_t> sort_indexes(std::vector<T> const &v)
{
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  
  return idx;
}

std::vector<double> getRanks(std::vector<double> const &absoluteValues)
{
  std::vector<double> ranks(absoluteValues.size());
  
  size_t i = 0;
  while (i < absoluteValues.size())
  {
    size_t j = i + 1;
    while (j < absoluteValues.size())
    {
      if(absoluteValues[i] != absoluteValues[j])
      {
        break;
      }
      j++;
    }
    for(size_t k = i; k <= j-1; k++)
    {   
      ranks[k] = 1 + (double)(i + j-1)/(double)2;
    }
    i = j;
  }
  return ranks;
}

std::vector<double> getUniq(std::vector<double> u)
{
  // using default comparison:
  std::vector<double>::iterator it;
  it = std::unique (u.begin(), u.end());
  u.resize( std::distance(u.begin(),it) );
  return u;
}

std::vector<double> getCounts(std::vector<double> const &valuesOrd)
{
  std::vector<double> u = getUniq(valuesOrd);
  std::vector<double> u_counts(u.size(),0);
  size_t k=0;
  double prev=valuesOrd[0];
  
  for(size_t i=0;i<valuesOrd.size();i++)
  {
    if(prev==valuesOrd[i])
    {
      u_counts[k]++;
    } else {
      k++;
      u_counts[k]++;
      prev=valuesOrd[i];
    }
  }
  return u_counts;
}

double getSigma(std::vector<double>& u_counts,double n1, double n2)
{
  double nties = 0;
  //int n = n1 + n2;
  
  if(u_counts.size() < n1 + n2)
  {
    for(std::vector<double>::iterator it = u_counts.begin(); it != u_counts.end(); ++it)
      nties += ( ((*it) * (*it) * (*it)) - *it );
  }
  
  return sqrt( (n1 * n2 / 12) * ( (n1 + n2 + 1) - nties / ((n1 + n2) * (n1 + n2 - 1))) );
}

double getPvalue(double& z)
{
  double p =1;
  if(z<0){
    p = gsl_cdf_gaussian_P(z,1)*2;
  } else {
    p = gsl_cdf_gaussian_Q(z,1)*2;
  }
  return p;
}

Rcpp::NumericVector subset(Rcpp::NumericVector const &v, Rcpp::NumericVector idx)
{
  Rcpp::NumericVector res(idx.size());
  size_t k=0;
  for(Rcpp::NumericVector::iterator it = idx.begin(); it != idx.end(); ++it, k++)
    res(k) = v[*it - 1];
  return res;
}

double avg(Rcpp::NumericVector const &v)
{
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double avg(std::vector<double> const &v)
{
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}


double MWUtest(Rcpp::NumericVector const &v1, Rcpp::NumericVector const &v2)
{
  double pval = 1;
  std::vector<double> values(v1.length());
  std::vector<double> valuesOrd(v1.length()+v2.length());
  std::vector<double> valuesRnk(v1.length()+v2.length());
  
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
    
    double mu = v1.size()*v2.size()/2, sig=0, z = 0;
    sig = getSigma(nties,v1.size(),v2.size());
    
    z = U1<U2 ? U1-mu : U2-mu;
    z = z<0 ? z+0.5 : z-0.5; // Continuity Correction
    z = z/sig;
    
    pval = getPvalue(z);
  }
  return pval;
}
