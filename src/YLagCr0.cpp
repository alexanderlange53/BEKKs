#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
arma::mat YLagCr0(arma::mat& r2, int& Tob, int& q, double& m_r2) {
  arma::mat help_m = arma::ones(r2.n_rows + q, q + 1);

  for (int j = 1; j <= q; j++) {
    int k = 2;
    for (int i = 0; i <= (q - 1); i++) {
      if (i < j) {
        help_m(i, j) = m_r2;
      } else {
        help_m(i, j) = r2(k);
        k+= 1;
      }
    }
  }

  return help_m.submat(0, 0, (q - 1), q);
}
