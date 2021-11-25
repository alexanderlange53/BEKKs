#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int indicatorFunction(arma::mat r,arma::mat signs){
  r=r.t();

  int indicator=1;
  int n=r.n_rows;
  for (int i=0; i<n; i++){
    if(signs(i,0)*r(i,0)<0){
      indicator=0;
    }
  }
  return indicator;
}
