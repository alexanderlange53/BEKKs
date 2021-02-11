#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double loglike_bekk(arma::vec theta, arma::mat r) {
// Log-Likelihood function

// convert to matrices
  int n = r.n_cols;
// Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
   }

   arma::mat A = arma::reshape(theta.subvec(index, (index + n^2) - 1 ), n, n);
   arma::mat G = arma::reshape(theta.subvec((index + n^2), numb_of_vars-1), n, n);

    Rcpp::Function f("valid_bekk");

// check constraints
    if (f(C, A, G) == 0) {
      return -1e25;
    }

// compute inital H
    arma::mat H = (r.t() * r) / r.n_rows;

    arma::mat CC  = C * C.t();
    arma::mat At  = A.t();
    arma::mat Gt  = G.t();

    double llv = -0.5 * arma::as_scalar(n * log(2 * M_PI) + log(arma::det(H)) + r.row(0) * arma::inv(H) * r.row(0).t());
    for (int i = 1; i < NoOBs; i++) {
      H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
      llv += -0.5 * arma::as_scalar(n * log(2 * M_PI) + log(arma::det(H)) + r.row(i) * arma::inv(H) * r.row(i).t());
    }
    return llv;
}

