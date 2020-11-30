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

// BHHH algorithm based on first order derivatives and outer product for Gaussian GARCH(p, q) model
// [[Rcpp::export]]
arma::vec BhhhGarch(arma::mat r2, int q, int p, arma::mat theta, arma::mat epsilon2, arma::mat Z, int Tob, int max_iter, double crit, double ucvar){

  // Args:
  //    - r2: squared returns
  //    - q: ARCH order
  //    - p: GARCH order
  //    - theta: parameter to optimize
  //    - epsilon2: GARCH process of squared returns
  //    - Z: Regressor matrix
  //    - Tob: Number of observations
  //    - max_iter: maximum nimber of algorithm iterations
  //    - crit: accuracy of algorithm
  //    - ucvar: unconditional variance

  // Generating candidate step width
  arma::vec step_width = arma::linspace(1, 0, 100);
  int quit = 0;
  int iter = 1;

  while (quit == 0 && iter < max_iter) {

    arma::mat score_function = ScoreGarch(epsilon2, Z, Tob, q, p, theta, ucvar);

    arma::mat out_prod_grad = score_function * score_function.t();

    arma::mat sum_grad = arma::sum(score_function, 1);

    arma::mat theta_hat = arma::zeros(theta.n_rows, 100);
    for (int i = 0; i < 100; i++) {
      theta_hat.col(i) = theta + step_width(i) * out_prod_grad.i() * sum_grad;
    }

    theta_hat = arma::abs(theta_hat);

    arma::vec candidates = arma::zeros(100);

    candidates(0) = LikelihoodGarch(Z, Tob, q, p, theta, epsilon2, ucvar);

    int quit_inner = 0;
    int iter_inner = 1;

    while (iter_inner < 100 && quit_inner == 0) {
      candidates(iter_inner) = LikelihoodGarch(Z, Tob, q, p, theta_hat.col(iter_inner), epsilon2, ucvar);
      iter_inner += 1;
    }

    int min_lik = candidates.index_min();

    double best_candidate = candidates.min();

    if (pow(best_candidate - candidates(0), 2) / abs(candidates(0)) <  crit) {
      quit = 1;
    } else {
      theta = theta_hat.col(min_lik);
      iter += 1;
    }
  }

  double lik_optim = LikelihoodGarch(Z, Tob, q, p, theta, epsilon2, ucvar);

  arma::vec theta_out = arma::vectorise(theta);
  int sz = theta_out.size();
  theta_out.resize(sz+1);
  theta_out(sz) = lik_optim;

  return theta_out;
}
