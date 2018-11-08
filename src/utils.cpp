// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <functional>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

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

typedef std::vector<uword> index_vec;
typedef std::pair<index_vec, std::vector<uword>> setPair;


struct Counter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type&) { ++count; }
  size_t count = 0;
};

uword count_overlap(index_vec V1, index_vec V2) {
  //V1 and V2 are sorted vectors
  Counter c;
  std::set_intersection(V1.begin(), V1.end(),
                        V2.begin(), V2.end(),
                        std::back_inserter(c));
  return c.count;
}

double jacc_pairs(setPair B1, setPair B2) {

  uword common_pairs = count_overlap(B1.first, B2.first) * 
                          count_overlap(B1.second, B2.second);
  
  uword total_pairs = B1.first.size() * B1.second.size() + 
                      B2.first.size() * B2.second.size() - 
                          common_pairs;
  
  
  return (double)common_pairs/(double)total_pairs;  
  
}

// [[Rcpp::export]]
NumericMatrix jacc_matrix_c(List bimod_list) {
  
  uword n = bimod_list.size();
  std::vector<setPair> bimods(n);
  
  for(uword i=0; i < n; i++) {
    bimods[i] = setPair(as<index_vec>(as<List>(bimod_list[i])["x"]),
                        as<index_vec>(as<List>(bimod_list[i])["y"]));
    std::sort(bimods[i].first.begin(), bimods[i].first.end());
    std::sort(bimods[i].second.begin(), bimods[i].second.end());
  }
  
  NumericMatrix res(n, n);
  
  for(uword i=0; i < n; i++) {
    res(i, i) = 1;
    for(uword j=0; j < i; j++) {
      res(i, j) = res(j, i) = jacc_pairs(bimods[i], bimods[j]);
    }
  }
  
  return res;
}

// [[Rcpp::export]]
double effective_num_c(List bimod_list) {
  uword n = bimod_list.size();
  std::vector<setPair> bimods(n);
  
  for(uword i=0; i < n; i++) {
    bimods[i] = setPair(as<index_vec>(as<List>(bimod_list[i])["x"]),
                        as<index_vec>(as<List>(bimod_list[i])["y"]));
  }
  
  std::unordered_map<uword, index_vec> ind_by_x, ind_by_y;

  for(uword i=0; i < n; i++) {
    for(uword x : bimods[i].first) {
      ind_by_x[x].push_back(i);
    }
    
    for(uword y : bimods[i].second) {
      ind_by_y[y].push_back(i);
    }
  }
  
  double total_overlap, expected_overlap;

  total_overlap = 0;
  
  for(const setPair& b: bimods) {
    
    expected_overlap = 0;
      
    for(uword x: b.first) {
      for(uword y: b.second) {
         double counts = count_overlap(ind_by_x.at(x), ind_by_y.at(y));
         expected_overlap += 1/counts;
      }
    }
  
  expected_overlap /= (b.first.size() * b.second.size());
  total_overlap += expected_overlap;
  }
  
  return total_overlap;
}