#include <Rcpp.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void updateColumnsInPlace(NumericMatrix &A,
                          const IntegerVector &ind, 
                         const NumericMatrix &B) {
  for(int i=0; i < ind.length(); i++) {
    A.column(ind[i]-1) = B.column(i);
  }
}
