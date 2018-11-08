#include <unordered_map>
#include <vector>
#include <RcppArmadillo.h>
#include <functional>
#include "lru_cache.cpp"
// we only include RcppArmadillo.h which pulls Rcpp.h in for us


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

const uword NUMERIC_SIZE = sizeof(double);
  
mat scale(const mat& X) {
  rowvec mu = mean(X);
  rowvec sigma = stddev(X);
  return (X.each_row() - mu).each_row()/sigma;
}

class CorBackend {

public:
  mat X, Y;
  uword dx, dy, n;

  typedef cache::lru_cache_vec<uword, double> cache_t;
  //The cache limit
  cache_t cache_x, cache_y;
  
  CorBackend(const mat &X, const mat &Y,
             uword cache_size_x, uword cache_size_y) :
    X(scale(X)), Y(scale(Y)),
    dx(X.n_cols), dy(Y.n_cols), n(Y.n_rows),
    cache_x(cache_size_x/(NUMERIC_SIZE*dy), dy), cache_y(cache_size_y/(NUMERIC_SIZE*dx), dx) {
      
    }
  
  mat calc_cor_x(uvec A) {
    return Y.t() * X.cols(A)/(n-1);
  }
  
  mat calc_cor_y(uvec A) {
    return X.t() * Y.cols(A)/(n-1);
  }
  
  mat calc_cor_memoized(const uvec& A, 
                 uword val_len,
                 cache_t& cache,
                 const std::function<mat(uvec)>& calc_cor) {
    uword m = A.n_elem;
    mat R(val_len, m);
    std::vector<uword> ind_new, A_new;
    
    for(uword i=0; i < m; i++) {
      uword a = A(i);
      if(cache.exists(a)) {
        R.col(i) = cache.get(a);
      } else {
        ind_new.push_back(i);
        A_new.push_back(a);
      }
    }
    
    //Compute the new correlations
    R.cols(conv_to<uvec>::from(ind_new)) = calc_cor(conv_to<uvec>::from(A_new));
    
    for(uword i=0; i < ind_new.size(); i++) {
      cache.put(A_new[i], R.col(ind_new[i]));
    }
    
    return R;
  }
  
  mat get_cor_x(const uvec& A) {
    return calc_cor_memoized(A, dy, cache_x, [this](uvec A) -> mat {
      return this->calc_cor_x(A);
    });
  }
  
  mat get_cor_y(const uvec& A) {
    return calc_cor_memoized(A, dx, cache_y, [this](uvec A) -> mat {
      return this->calc_cor_y(A);
    });
  }
};


mat get_cor_export(CorBackend *bk, uvec A) {
  uword dx = bk->dx;
  if(A.min() > dx) {
    A.transform([dx](uword i) { return (i - dx - 1);});
    return bk->get_cor_y(A);
  } else {
    A.transform([](uword i) { return (i - 1);});
    return bk->get_cor_x(A); 
  }
} 

vec get_sq_tstat(CorBackend *bk, const uvec& A) {
 return sum(square(get_cor_export(bk, A)), 1);
} 

mat get_scale_mat_export(CorBackend *bk, uvec A) {
  uword dx = bk->dx;
  if(A.min() > dx) {
    A.transform([dx](uword i) { return (i - dx - 1);});
    return bk->Y.cols(A);
  } else {
    A.transform([](uword i) { return (i - 1);});
    return bk->X.cols(A); 
  }
}

RCPP_MODULE(cbase) {
  using namespace Rcpp;

  class_<CorBackend>("CorBackend")
//     // expose the default constructor
      .constructor<arma::mat, arma::mat, arma::uword,  arma::uword>()
      .method("getCor", &get_cor_export , "Get the correrlations")
      .method("getSqTstat", &get_sq_tstat , "Get the sum of squared test statistic")
      .method("getScaledMat", &get_scale_mat_export, "Get the scaled matrices")
      .field("dx", &CorBackend::dx)
      .field("dy", &CorBackend::dy)
      .field("n", &CorBackend::n)
//     .method("set", &CorBackend::set     , "set the message")
     ;
}

