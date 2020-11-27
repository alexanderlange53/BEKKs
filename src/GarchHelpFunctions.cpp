#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Generates lag matrix (regressor matrix) of recursive volatilities
// [[Rcpp::export]]
arma::mat SigmaLagCr(arma::mat y, int Tob, int q, int p, double ucvar){

  // Args:
  //    - y: vector of conditional GARCH variances
  //    - Tob: Number of observations
  //    - q: ARCH order
  //    - p: GARCH order
  //    - ucvar: unconditional variance

  arma::mat SigmaLag = arma::zeros(Tob-q,p);

  int k = 0;

  for (int j = 0; j < p; ++j) {
    k = 0;
    for (int i = 0; i < (Tob-q); ++i) {
      if (j >= i) {
        SigmaLag(i, j) = ucvar;
      } else{
        SigmaLag(i, j) = y(0, k);
        k += 1;
      }
    }
  }

  return SigmaLag;
}

// Generates vector of GARCH variances (Sigma_e)
// [[Rcpp::export]]
arma::mat GarchVariance(int Tob, int q, int p, double ucvar, arma::mat theta, arma::mat Z){

  // Args:
  //    - Tob: Number of observations
  //    - q: ARCH order
  //    - p: GARCH order
  //    - ucvar: unconditional variance
  //    - theta: parameter to optimize
  //    - Z: Regressor matrix

  arma::mat vvec = arma::zeros(1, Tob-q);

  arma::mat e = ucvar * arma::ones(1, 1, p);
  vvec = arma::join_horiz(e, vvec);

  int g = theta.n_rows;

  for (int i = p; i < (Tob-q+p); i++) {
    arma::mat b = arma::reverse(vvec.cols((i-p), (i-1)), 1);
    b = b * theta.rows(q+1, g-1);
    vvec.col(i) = Z.row(i-p) * theta.rows(0, q) + b;
  }

  return vvec.cols(p, (Tob-q+p-1));
}

// Score function of Gaussian log-likelihood for GARCH(p, q) model
// [[Rcpp::export]]
arma::mat ScoreGarch(arma::mat epsilon2, arma::mat Z, int Tob, int q, int p, arma::mat theta, double ucvar) {

  // Args:
  //    - epsilon2: GARCH process of squared returns
  //    - Z: Regressor matrix
  //    - Tob: Number of observations
  //    - q: ARCH order
  //    - p: GARCH order
  //    - theta: parameter to optimize
  //    - ucvar: unconditional variance

  // GARCH variance sigma^2
  arma::mat vvec = GarchVariance(Tob, q, p, ucvar, theta, Z);

  // Lag matrix for recursive  volatility
  arma:: mat sigma_lag_mat = SigmaLagCr(vvec, Tob, q, p, ucvar);

  double aa;
  aa = arma::as_scalar(ucvar/(1 - arma::sum(theta.rows((q + 1), (theta.n_rows-1)))));
  double aa2;
  aa2 = arma::as_scalar(1/(1 - arma::sum(theta.rows((q + 1), (theta.n_rows-1)))));

  arma::mat gz = join_horiz(Z, sigma_lag_mat);

  arma::mat full_regressor_mat = arma::join_horiz(Z, sigma_lag_mat);
  int nparam = theta.n_rows;

  arma::mat vabl = arma::zeros(Tob-q, nparam);

  arma::mat e = aa * arma::ones(p, nparam);
  e.col(0) = aa2 * arma::ones(p, 1);
  vabl = join_vert(e, vabl);

  for (int i = p; i < (Tob-q+p); i++) {
    arma::mat b = vabl.rows((i - p), (i-1));
    b = arma::repmat(theta.rows(q+1, nparam-1), 1, b.n_cols) % b;
    if (p > 1) {
      b = arma::sum(b);
    }

    vabl.row(i) = gz.row(i-p) + b;
  }

  vabl = vabl.rows(p, (Tob-q+p-1));

  arma::mat score1 = vabl/arma::repmat(vvec.t(), 1, vabl.n_cols);
  arma::mat score2 = arma::pow(vvec, 2);
  score2 = epsilon2/score2;
  score2 = arma::repmat(score2.t(), 1, vabl.n_cols) % vabl;

  arma::mat ltv = 0.5 * (score2 - score1);

  return ltv.t();
}

// Gaussian log-likelihood for GARCH(p, q) model
// [[Rcpp::export]]
double LikelihoodGarch(arma::mat Z, int Tob, int q, int p, arma::mat theta, arma::mat epsilon2, double ucvar) {

  // Args:
  //    - Z: Regressor matrix
  //    - Tob: Number of observations
  //    - q: ARCH order
  //    - p: GARCH order
  //    - epsilon2: GARCH process of squared returns
  //    - ucvar: unconditional variance

  arma::mat sigma2 = GarchVariance(Tob, q, p, ucvar, theta, Z);

  arma::mat term2 = epsilon2 / sigma2;

  double loglik = 0.5 * (Tob - q) * log(2 * M_PI) + arma::as_scalar(0.5 * arma::sum(log(sigma2), 1)) + 0.5 * arma::as_scalar(arma::sum(term2, 1));
  return loglik;
}
