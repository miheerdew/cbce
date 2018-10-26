// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do

//' Update the colmns ind of matrix A in place by matrix B
//' 
//' @export
// [[Rcpp::export]]
void updateColumnsInPlace(arma::mat &A,
                          arma::uvec ind, 
                          const arma::mat &B) {
  ind.transform([](arma::uword i) -> arma::uword { return (i - 1); });
  A.cols(ind) = B;
}
