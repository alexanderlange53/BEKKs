#include <RcppArmadillo.h>



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


#ifndef Indicator_Functions_H
#define Indicator_Functions_H

// [[Rcpp::export]]
inline int indicatorFunction(arma::mat r, arma::mat signs){
  r = r.t();

  int indicator = 1;
  int n = r.n_rows;
  for (int i = 0; i<n; i++){
    if(arma::as_scalar(signs.row(i)) * arma::as_scalar(r.row(i)) < 0){
      indicator = 0;
    }
  }
  return indicator;
}

#endif
