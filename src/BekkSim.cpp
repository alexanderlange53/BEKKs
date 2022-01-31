#include <RcppArmadillo.h>

 #include "IndicatorFunctions.h"

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

// [[Rcpp::export]]
arma::mat simulate_bekka_c(arma::vec theta, const int NoObs, const int n, arma::vec signs, double expected_signs){

  // Length of each series
  arma::mat series = arma::zeros(n, NoObs);
  int numb_of_vars = 3 * pow(n, 2) + n*(n + 1)/2;
  arma::mat C = arma::zeros(n, n);
  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }


  arma::mat innovations(n, NoObs, arma::fill::randn);


  // double exp_indicator_value = 0;
  // for (int i=0; i<NoObs;i++){
  //   exp_indicator_value+=indicatorFunction(innovations.row(i),signs);
  // }
  //exp_indicator_value=exp_indicator_value/NoObs;

  arma::mat C_full = C * C.t();
  arma::mat A = arma::reshape(theta.subvec(index, (index + pow(n, 2) - 1)), n, n);
  arma::mat At = A.t();
  arma::mat B = arma::reshape(theta.subvec((index + pow(n, 2)), index + 2*pow(n, 2) - 1), n, n);
  arma::mat Bt= B.t();
  arma::mat G = arma::reshape(theta.subvec((index + 2*pow(n, 2)), numb_of_vars - 1), n, n);
  arma::mat Gt= G.t();


  // unconditional variance replaced pow(n,2) by n
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(n, 2), pow(n, 2)) - arma::trans(arma::kron(A, A)) - expected_signs*arma::trans(arma::kron(B, B))- arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), n, n);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();


  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
    H = C_full + At * series.col(i-1) * series.col(i-1).t() * A + indicatorFunction(series.col(i-1).t(),signs)*Bt * series.col(i-1) * series.col(i-1).t() * B + Gt * H * G;

    h_dec = arma::chol(H).t();
    series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}

// [[Rcpp::export]]
arma::mat simulate_dbekk_c(arma::vec theta, const int NoObs, const int N){

  // Length of each series
  arma::mat series = arma::zeros(N, NoObs);
  int numb_of_vars = 2 * pow(N, 2) + N*(N + 1)/2;
  arma::mat C = arma::zeros(N, N);
  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }

  arma::mat C_full = C * C.t();
  arma::mat A = arma::diagmat(theta.rows(N * (N+1)/2, N + (N * (N + 1)/2) - 1));
  arma::mat G = arma::diagmat(theta.rows(N + N * (N + 1)/2, 2*N + (N * (N + 1)/2) - 1));


  // unconditional variance
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(N, 2), pow(N, 2)) - arma::trans(arma::kron(A, A)) - arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), N, N);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();

  arma::mat innovations(N, NoObs, arma::fill::randn);
  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
    H = C_full + A * series.col(i-1) * series.col(i-1).t() * A + G * H * G;
    h_dec = arma::chol(H).t();
    series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}

// [[Rcpp::export]]
arma::mat simulate_dbekka_c(arma::vec theta, const int NoObs, const int N, arma::vec signs, double expected_signs){

  // Length of each series
  arma::mat series = arma::zeros(N, NoObs);
  int numb_of_vars = 3 * pow(N, 2) + N*(N + 1)/2;
  arma::mat C = arma::zeros(N,N);
  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }


  arma::mat innovations(N, NoObs, arma::fill::randn);


  // double exp_indicator_value = 0;
  // for (int i=0; i<NoObs;i++){
  //   exp_indicator_value+=indicatorFunction(innovations.row(i),signs);
  // }
  //exp_indicator_value=exp_indicator_value/NoObs;

  arma::mat C_full = C * C.t();
  arma::mat A = arma::diagmat(theta.rows(N * (N+1)/2, N + (N * (N + 1)/2) - 1));
  arma::mat B = arma::diagmat(theta.rows(N + N * (N + 1)/2, 2*N + (N * (N + 1)/2) - 1));
  arma::mat G = arma::diagmat(theta.rows(2*N + N * (N + 1)/2, 3*N + (N * (N + 1)/2) - 1));

  // unconditional variance replaced pow(n,2) by n
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(N, 2), pow(N, 2)) - arma::trans(arma::kron(A, A)) - expected_signs*arma::trans(arma::kron(B, B))- arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), N, N);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();


  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
    H = C_full + A * series.col(i-1) * series.col(i-1).t() * A + indicatorFunction(series.col(i-1).t(),signs)*B * series.col(i-1) * series.col(i-1).t() * B + G * H * G;

    h_dec = arma::chol(H).t();
    series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}




// [[Rcpp::export]]
arma::mat simulate_sbekk_c(arma::vec theta, const int NoObs, const int N){

  // Length of each series


  double a = theta[0];
  double b = theta[1];

  // unconditional variance
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(N, 2), pow(N, 2)) - arma::trans(arma::kron(A, A)) - arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), N, N);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();

  arma::mat innovations(N, NoObs, arma::fill::randn);
  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
    H = C_full + A * series.col(i-1) * series.col(i-1).t() * A + G * H * G;
    h_dec = arma::chol(H).t();
    series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}

// [[Rcpp::export]]
arma::mat simulate_sbekka_c(arma::vec theta, const int NoObs, const int N, arma::vec signs, double expected_signs){

  // Length of each series
  arma::mat series = arma::zeros(N, NoObs);
  int numb_of_vars = 3 * pow(N, 2) + N*(N + 1)/2;
  arma::mat C = arma::zeros(N,N);
  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }


  arma::mat innovations(N, NoObs, arma::fill::randn);


  // double exp_indicator_value = 0;
  // for (int i=0; i<NoObs;i++){
  //   exp_indicator_value+=indicatorFunction(innovations.row(i),signs);
  // }
  //exp_indicator_value=exp_indicator_value/NoObs;

  arma::mat C_full = C * C.t();
  arma::mat A = arma::diagmat(theta.rows(N * (N+1)/2, N + (N * (N + 1)/2) - 1));
  arma::mat B = arma::diagmat(theta.rows(N + N * (N + 1)/2, 2*N + (N * (N + 1)/2) - 1));
  arma::mat G = arma::diagmat(theta.rows(2*N + N * (N + 1)/2, 3*N + (N * (N + 1)/2) - 1));

  // unconditional variance replaced pow(n,2) by n
  arma::mat Uncond_var = arma::reshape(arma::inv(arma::eye(pow(N, 2), pow(N, 2)) - arma::trans(arma::kron(A, A)) - expected_signs*arma::trans(arma::kron(B, B))- arma::trans(arma::kron(G, G))) * arma::vectorise(C_full), N, N);

  arma::mat H = Uncond_var;
  arma::mat h_dec = arma::chol(H).t();


  series.col(0) = h_dec * innovations.col(0);

  for(int i = 1; i < NoObs; i++){
    H = C_full + A * series.col(i-1) * series.col(i-1).t() * A + indicatorFunction(series.col(i-1).t(),signs)*B * series.col(i-1) * series.col(i-1).t() * B + G * H * G;

    h_dec = arma::chol(H).t();
    series.col(i) = h_dec * innovations.col(i);
  }

  return series.t();
}

