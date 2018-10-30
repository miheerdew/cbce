#include <unordered_map>
#include <vector>
#include <RcppArmadillo.h>
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

typedef std::unordered_map<uword, uword> cached_db;
typedef std::unordered_map<uword, int> timestamp_db;
  
mat get_cors(const uvec& A,
             mat& cache,
             cached_db& cached,
             timestamp_db& tstamp,
             int& elapsed,
             const std::function<mat(uvec)>& calc_cor) {
    // Calculate the corelations for A.
    // A p x m matrix  where m = len(A)
    uword m = A.n_elem;
    uword p = cache.n_rows;
    
    std::vector<uword> A_new(m), A_old(m), old_pos(m); 
    
    std::for_each(A.begin(), A.end(), [&](uword a) {
      auto search = cached.find(a);
      if(search != cached.end()) {
        A_old.push_back(a);
        old_pos.push_back(search->second);
      } else {
        A_new.push_back(a);
      }
    });
      
    elapsed += 1;
    
    mat R(p, m);
    
    if (A_old.size() > 0 ) {
      //Assign the correlations already in cache.
      R.cols(conv_to< uvec >::from(A_old)) = cache.cols(conv_to< uvec >::from(old_pos));
      std::for_each(A_old.begin(), A_old.end(), [&](uword a) {
        //Assign the latest time stamp
        tstamp[a] = elapsed;
      });
    }
        
    if (A_new.size() > 0) {
          uvec A_new_v = conv_to< uvec >::from(A_new);
          // Calculate the new correlations
          R.cols(A_new_v) = calc_cor(A_new);
          
          // Find the indices which must be removed 
          //Populate the cache with these correlations.
    }
    
    return R;
}

class CorBackend {

public:
  mat X, Y;
  uword dx, dy, n;

  //The cache limit
  uword cache_lim_x, cache_lim_y;
  mat cache_x, cache_y;

  //What indices are cached:
  //cached_db cached_x, cached_y;

  //Time stamps
  //timestamp_db tstamp_x, tstamp_y;

  CorBackend(const mat &X, const mat &Y,
             uword cache_size_x, uword cache_size_y) :
    X(scale(X)), Y(scale(Y)),
    dx(X.n_cols), dy(Y.n_cols), n(Y.n_rows),
    cache_lim_x(cache_size_x/(NUMERIC_SIZE*dy)), cache_lim_y(cache_size_y/(NUMERIC_SIZE)*dx),
    cache_x(dy, cache_lim_y), cache_y(dx, cache_lim_y) {}
  
  mat calc_cor_x(uvec A) {
    return Y.t() * X.cols(A)/(n-1);
  }
  
  mat calc_cor_y(uvec A) {
    return X.t() * Y.cols(A)/(n-1);
  }
  
  mat get_cor_x(uvec A) {
    return calc_cor_x(A);
  }
  
  mat get_cor_y(uvec A) {
    return calc_cor_y(A);
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

RCPP_MODULE(cbase) {
  using namespace Rcpp;

  class_<CorBackend>("CorBackend")
//     // expose the default constructor
      .constructor<arma::mat, arma::mat, arma::uword,  arma::uword>()
      .method("getCor", &get_cor_export , "Get the correrlations")
//     .method("set", &CorBackend::set     , "set the message")
     ;
}

