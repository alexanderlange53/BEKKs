#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::mat simulate_bekk_c(arma::vec theta, const int NoObs, const int n){

  // Length of each series
  arma::mat series = arma::zeros(n, NoObs);
  int numb_of_vars = 2 * pow(n, 2) + n*(n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }

  arma::mat C_full = C * C.t();
  arma::mat A = arma::reshape(theta.subvec(index, (index + pow(n, 2) - 1)), n, n);
  arma::mat At = A.t();
  arma::mat G = arma::reshape(theta.subvec((index + pow(n, 2)), numb_of_vars - 1), n, n);
  arma::mat Gt= G.t();


  // unconditional variance
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(n, 2), pow(n, 2)) - arma::trans(arma::kron(A, A)) - arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), n, n);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();

  arma::mat innovations(n, NoObs, arma::fill::randn);
  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
      H = C_full + At * series.col(i-1) * series.col(i-1).t() * A + Gt * H * G;
      h_dec = arma::chol(H).t();
      series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}
