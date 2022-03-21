#include <RcppArmadillo.h>

#include "IndicatorFunctions.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::mat indicatorFunction(arma::mat r, arma::mat signs){
  r = r.t();
  arma::mat indicator = r;
  //int indicator = 1;
  int n = r.n_rows;
  for (int i = 0; i<n; i++){
    if(arma::as_scalar(signs.row(i)) * arma::as_scalar(r.row(i)) < 0){
      indicator[i] = 0;
    }
  }
  return indicator;
}

// [[Rcpp::export]]
arma::mat expected_indicator_value(arma::mat r, arma::mat signs){
  int N =r.n_rows;
  int n =r.n_cols;
  arma::mat S2 = r.t() * r / N;

  arma::mat S1 = arma::zeros(n,n);
  for (int i=0; i<N;i++){
    S1+=indicatorFunction(r.row(i),signs)*indicatorFunction(r.row(i),signs).t();
  }
  arma::mat exp_indicator_value=(S1/N)/S2;
  return exp_indicator_value;
}

// [[Rcpp::export]]
arma::mat elimination_mat(const int& n) {
  // Generates an elimination matrix for size 'n'
  int n1 = n * (n + 1) / 2;
  int n2 = pow(n, 2);

  arma::mat init = arma::eye(n1, n1);
  int oes = 1;

  arma::mat eli  = init.col(0);
  int block = n;

  while (eli.n_cols < n2) {
    if (eli.n_cols == 1) {
      eli = init.cols(0, block-1);
    } else {
      eli = arma::join_horiz(eli, init.cols(0, block-1));
    }

    if (init.n_cols > 1) {
      init = init.cols(block, init.n_cols-1);
    }

    eli = arma::join_horiz(eli, arma::zeros(eli.n_rows, oes));

    oes += 1;

    block -= 1;
  }

  return eli.cols(0, eli.n_cols-n-1);
}

// [[Rcpp::export]]
arma::mat commutation_mat(const int& n) {
  // generates a (square) commutation matrix for 'n'
  arma::mat K = arma::zeros(pow(n, 2), pow(n, 2));

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      K(i + n*(j - 1)-1, j + n*(i - 1)-1) = 1;
    }
  }
  return K;
}

// [[Rcpp::export]]
arma::mat duplication_mat(const int& n) {
  // Generates a duplication matrix for size 'n'
  int n2 = pow(n, 2);

  arma::mat el = elimination_mat(n);
  arma::mat co = commutation_mat(n);
  arma::mat m = arma::eye(n2, n2) + co;

  arma::mat dup = m*el.t()*arma::inv(el*m*el.t());

  return dup;
}

// [[Rcpp::export]]
arma::mat inv_gen(const arma::mat& m) {
  // Checks if a matrix is positive definit and calculates
  // the inverse or generalized inverse

  if (m.is_sympd() == TRUE) {
    return m.i();
  } else {
    return arma::pinv(m);
  }
}
