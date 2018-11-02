// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do

//' Update the colmns ind of matrix A in place by matrix B
//' 
// [[Rcpp::export]]
void updateColumnsInPlace(arma::mat &A,
                          arma::uvec ind, 
                          const arma::mat &B) {
  ind.transform([](arma::uword i) -> arma::uword { return (i - 1); });
  A.cols(ind) = B;
}

typedef std::vector<double> std_vec;

std::vector<uword> fast_bh(const std_vec& Tstat, 
                           double alpha, 
                           const std::function<double(double)>& quantile,
                           const std::function<bool(double, double)>& cmp) {
  uword m = Tstat.size();
  
  std_vec tmp(Tstat);
  std_vec::iterator it_start, it_end;
  
  it_start = tmp.begin();
  it_end = tmp.end();
  
  uword d, d_old;
  double quant;
  
  d = std::distance(it_start, it_end);

  do {
    d_old = d;
    quant = quantile(alpha*d/m);
    it_end = std::partition(it_start, it_end, 
                            [quant, cmp](double t) { return cmp(t, quant); });
    d = std::distance(it_start, it_end);
  } while (d != d_old);
  
  std::vector<uword> indices;
  
  if(d > 0) {
    for(uword i=0; i < m; i++) {
      if(cmp(Tstat[i], quant)) {
        indices.push_back(i);
      }
    }
  }
  
  return indices;
}


//'
//'@export
// [[Rcpp::export]]
IntegerVector fast_bh_beta(NumericVector Tstat, double alpha, 
                           double shape1, double shape2,
                           bool lower=false) {
  std::vector<uword> res;
  
  if(lower) {
    res = fast_bh(as<std_vec>(Tstat), 
                        alpha,
                        [shape1, shape2, lower](double p) {
                          return R::qbeta(p, shape1, shape2, lower, false);
                        }, 
                        std::less_equal<double>());
  } else {
    res = fast_bh(as<std_vec>(Tstat), 
                  alpha,
                  [shape1, shape2, lower](double p) {
                    return R::qbeta(p, shape1, shape2, lower, false);
                  }, 
                  std::greater_equal<double>());
  }
  IntegerVector resR = wrap(res);
  return (resR + 1);
}
  
//'
//'@export
// [[Rcpp::export]]
IntegerVector fast_bh_chisq(NumericVector Tstat, 
                             double alpha, 
                             double df,
                             bool lower=false) {
    std::vector<uword> res;
    if(lower) {
      res = fast_bh(as<std_vec>(Tstat), 
                  alpha,
                  [df, lower](double p) {
                    return R::qchisq(p, df, lower, false);
                  },
                  std::less_equal<double>());
    } else {
      res = fast_bh(as<std_vec>(Tstat), 
                  alpha,
                  [df, lower](double p) {
                    return R::qchisq(p, df, lower, false);
                  },
                  std::greater_equal<double>());
   }
  
   IntegerVector resR = wrap(res);
   return (resR + 1);
}

