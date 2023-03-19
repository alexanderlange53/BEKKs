#include <RcppArmadillo.h>
#include "IndicatorFunctions.h"




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]



// [[Rcpp::export]]
bool valid_sbekk(arma::mat C, double a, double g){
  double sum=a+g;
  if(sum >= 1){
    return false;
  }


  int n =C.n_cols;



  for (int i=0; i<n;i++){
    if(C(i,i)<=0){
      return false;
    }
  }

  if(a<=0 || g<=0) {
    return false;
  }
  else{
    return true;
  }
}

// [[Rcpp::export]]
bool valid_asymm_sbekk(arma::mat C, double a, double b, double g, arma::mat r, arma::mat signs){
  double exp_indicator_value = expected_indicator_value(r,signs);
  double sum=a+g+b*exp_indicator_value;
  if(sum >= 1){
    return false;
  }
  int n =C.n_cols;
  for (int i=0; i<n;i++){
    if(C(i,i)<=0){
      return false;
    }
  }

  if(a<=0 || b<=0  || g<=0) {
    return false;
  }
  else{
    return true;
  }
}


// [[Rcpp::export]]
double loglike_sbekk(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int NofVarsC = theta.n_rows -2;
  double a = theta[NofVarsC];
  double g = theta[NofVarsC+1];
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat CC  = C * C.t();

  // check constraints
  if (valid_sbekk(C,a,g) == false) {
    return -1e25;
  }

  // compute inital H
  arma::mat  sample_covariance = (r.t() * r) / r.n_rows;
  arma::mat  H = sample_covariance;


  double llv = arma::as_scalar(log(arma::det(H)) + r.row(0) * arma::inv(H) * r.row(0).t());
  for (int i = 1; i < NoOBs; i++) {
    H = CC + a * r.row(i - 1).t() * r.row(i - 1)  + g * H ;
    llv += arma::as_scalar(log(arma::det(H)) + r.row(i) * arma::inv(H) * r.row(i).t());
  }


  return -0.5 * n * NoOBs * log(2 * M_PI) - 0.5 * llv;
}

// [[Rcpp::export]]
double loglike_asymm_sbekk(const arma::vec& theta, const arma::mat& r, arma::mat& signs) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;

  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat CC  = C * C.t();

  int NofVarsC = theta.n_rows -3;

  double a = theta[NofVarsC];
  double b = theta[NofVarsC+1];
  double g = theta[NofVarsC+2];

  double expected_indicator = expected_indicator_value(r, signs);

  // check constraints
  if (valid_asymm_sbekk(C,a,b,g,r,signs) == false) {
    return -1e25;
  }

  // compute inital H
  arma::mat  sample_covariance = (r.t() * r) / r.n_rows;
  arma::mat  H = sample_covariance;


  double llv = arma::as_scalar(log(arma::det(H)) + r.row(0) * arma::inv(H) * r.row(0).t());
  for (int i = 1; i < NoOBs; i++) {
    H = CC + (a+b*indicatorFunction(r.row(i - 1),signs)) * r.row(i - 1).t() * r.row(i - 1)  + g * H ;
    llv += arma::as_scalar(log(arma::det(H)) + r.row(i) * arma::inv(H) * r.row(i).t());
  }


  return -0.5 * n * NoOBs * log(2 * M_PI) - 0.5 * llv;
}


// [[Rcpp::export]]
arma::mat score_sbekk(const arma::mat& theta, arma::mat& r) {

  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);


  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);


  // Length of each series
  int NoOBs = r.n_rows;
  int NofVarsC = theta.n_rows -2;
  double a = theta[NofVarsC];
  double g = theta[NofVarsC+1];
  arma::mat C = arma::zeros(N, N);

  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat CC  = C * C.t();

  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model

  // Partial derivatives for initial period t = 1
  arma::mat sample_covar = r.t() * r / NoOBs;
  arma::mat ht = sample_covar;

  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t();

  arma::mat dHda = arma::zeros(N2, 1);

  arma::mat dHdg = arma::zeros(N2, 1);

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

  arma::mat ht_sqrt_inv = arma::inv(ht);
  for (int k = 0; k < theta.n_rows; k++) {

    arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
    //Rcpp::Rcout << ht_sqrt_inv;
    //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
    arma::mat mat_temp = dh * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dh * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, k) = -(0.5) * arma::sum(mat_temp.diag());
    //Rcpp::Rcout << gradients(0, k);
  }


  // Partial derivatives for period t >= 2

  for (int i = 1; i < r.n_rows; i++) {
    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t() + g * dHdc;



    dHda = arma::reshape( r.row(i-1).t() * r.row(i-1),N2,1) + g * dHda;

    //dHda =  arma::kron(arma::eye(N, N),at*r.row(i-1).t() *r.row(i-1)) + arma::kron(at * r.row(i-1).t() *r.row(i-1), arma::eye(N, N)) * commutation_mat(N)+ GGt * dHda;

    dHdg = arma::reshape(ht,N2,1) + g * dHdg;

    //dHdg = arma::kron(arma::eye(N, N),gt*ht) + arma::kron(gt * ht , arma::eye(N, N)) * commutation_mat(N) + GGt * dHdg;

    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

    ht = CC + a * r.row(i-1).t() * r.row(i-1) + g*ht;

    //ht_sqrt_inv = arma::inv(arma::real(arma::sqrtmat(ht)));
    arma::mat ht_sqrt_inv = inv_gen(ht);
    //et = ht_sqrt_inv * r.row(i).t();


    for (int k = 0; k < theta.n_rows; k++) {
      //arma::mat dh = arma::reshape(dHdtheta.row(k), N, N).t();
      arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
      //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
      arma::mat mat_temp = dh * ht_sqrt_inv - r.row(i).t() * r.row(i)* ht_sqrt_inv * dh * ht_sqrt_inv;
      gradients(i, k) = -(0.5) * arma::sum(mat_temp.diag());
    }
  }


  return gradients;
}

// [[Rcpp::export]]
arma::mat score_asymm_sbekk(const arma::mat& theta, arma::mat& r, arma::mat& signs) {

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);
  double expected_indicator = expected_indicator_value(r, signs);
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);


  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();




  // Length of each series
  int NoOBs = r.n_rows;
  int NofVarsC = theta.n_rows -3;
  double a = theta[NofVarsC];
  double b = theta[NofVarsC+1];
  double g = theta[NofVarsC+2];
  arma::mat C = arma::zeros(N, N);

  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat CC  = C * C.t();

  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model

  // Partial derivatives for initial period t = 1
  arma::mat sample_covar = r.t() * r / NoOBs;
  arma::mat ht = sample_covar;

  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t();

  arma::mat dHda = arma::zeros(N2, 1);
  arma::mat dHdb = arma::zeros(N2, 1);
  arma::mat dHdg = arma::zeros(N2, 1);

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdb, dHdg).t();

  arma::mat ht_sqrt_inv = arma::inv(ht);
  for (int k = 0; k < theta.n_rows; k++) {

    arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
    //Rcpp::Rcout << ht_sqrt_inv;
    //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
    arma::mat mat_temp = dh * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dh * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, k) = -(0.5) * arma::sum(mat_temp.diag());
    //Rcpp::Rcout << gradients(0, k);
  }


  // Partial derivatives for period t >= 2

  for (int i = 1; i < r.n_rows; i++) {
    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t() + g * dHdc;



    dHda = arma::reshape( r.row(i-1).t() * r.row(i-1),N2, 1) + g * dHda;
    dHdb = indicatorFunction(r.row(i - 1),signs)*arma::reshape( r.row(i-1).t() * r.row(i-1),N2, 1) + g * dHdb;

    dHdg = arma::reshape(ht,N2, 1) + g * dHdg;

    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdb, dHdg).t();

    ht = CC + (a+b*indicatorFunction(r.row(i - 1),signs)) * r.row(i - 1).t() * r.row(i - 1)  + g * ht ;


    //ht_sqrt_inv = arma::inv(arma::real(arma::sqrtmat(ht)));
    arma::mat ht_sqrt_inv = inv_gen(ht);
    //et = ht_sqrt_inv * r.row(i).t();


    for (int k = 0; k < theta.n_rows; k++) {
      //arma::mat dh = arma::reshape(dHdtheta.row(k), N, N).t();
      arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
      //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
      arma::mat mat_temp = dh * ht_sqrt_inv - r.row(i).t() * r.row(i)* ht_sqrt_inv * dh * ht_sqrt_inv;
      gradients(i, k) = -(0.5) * arma::sum(mat_temp.diag());
    }
  }

  return gradients;
}

// [[Rcpp::export]]
Rcpp::List  bhh_sbekk(arma::mat& r, const arma::mat& theta, int& max_iter, double& crit) {

  arma::vec steps = {9.9,9,8,7,6,5,4,3,2,1,0.5,0.25,0.1,0.01,0.005,
                     0.001,0.0005,0.0001,0.00005,0.00001,0};
  double step = 0.1;
  int count_loop = 0;
  arma::mat theta_candidate = theta;
  int exit_loop = 0;
  arma::vec lik_all(max_iter+1, arma::fill::zeros);
  lik_all(0) = loglike_sbekk(theta, r);


  while (count_loop < max_iter && exit_loop == 0) {
    arma::mat theta_loop = theta_candidate;
    arma::mat theta_temp = arma::zeros(theta_loop.n_rows, steps.n_elem);

    arma::mat score_function = score_sbekk(theta_loop, r);
    arma::mat outer_score = score_function.t() * score_function;

    arma::mat outer_score_inv = arma::inv(outer_score);
    arma::mat score_function_sum = arma::sum(score_function);

    double lik = loglike_sbekk(theta_loop, r);

    for (int i = 0; i < steps.n_elem; i++) {
      arma::vec temp = theta_candidate + step * steps(i) * outer_score_inv * score_function_sum.t();
      theta_temp.col(i) = temp;
    }


    arma::vec likelihood_candidates(steps.n_elem, arma::fill::zeros);
    likelihood_candidates(steps.n_elem - 1) = lik;

    int  j = steps.n_elem - 2;
    int exit_inner = 0;
    while (j >= 0 && exit_inner == 0) {
      likelihood_candidates(j) = loglike_sbekk(theta_temp.col(j), r);
      //if (likelihood_candidates(j+1) > likelihood_candidates(j)) {
      //  exit_inner = 1;
      //}
      j -= 1;
    }

    //return likelihood_candidates;

    //int max_index = arma::index_max(likelihood_candidates.subvec(j, (steps.n_elem -1))) + j;
    int max_index = arma::index_max(likelihood_candidates);
    //return max_index;
    double likelihood_best = likelihood_candidates(max_index);

    // exit criterion strange
    //if (pow(likelihood_best - likelihood_candidates(steps.n_elem -1), 2)/abs(likelihood_candidates(steps.n_elem -1)) < crit) {
    //  exit_loop = 1;
    //}
    if (likelihood_best < lik_all(count_loop)) {
      exit_loop = 1;
      count_loop += 1;
    } else if (pow(likelihood_best - likelihood_candidates(steps.n_elem -1), 2)/abs(likelihood_candidates(steps.n_elem -1)) < crit) {// if (pow(likelihood_best - lik_all(count_loop), 2)/abs(lik_all(count_loop)) < crit) {
      exit_loop = 1;
      count_loop += 1;
      theta_candidate = theta_temp.col(max_index);
      lik_all(count_loop) = likelihood_candidates(steps.n_elem -1);
    } else {
      theta_candidate = theta_temp.col(max_index);
      count_loop += 1;
      lik_all(count_loop) = likelihood_candidates(steps.n_elem -1);
    }
  }

  double likelihood_final = loglike_sbekk(theta_candidate, r);
  arma::mat score_final = score_sbekk(theta_candidate, r);
  arma::mat s1_temp = arma::diagmat(arma::inv(score_final.t() * score_final));
  arma::mat s1 = arma::sqrt(s1_temp.diag());

  arma::mat t_val = theta_candidate/s1;
  return Rcpp::List::create(Rcpp::Named("theta") = theta_candidate,
                            Rcpp::Named("t_val") = t_val,
                            Rcpp::Named("likelihood") = likelihood_final,
                            Rcpp::Named("iter") = count_loop,
                            Rcpp::Named("likelihood_iter") = lik_all);
}

// [[Rcpp::export]]
Rcpp::List  bhh_asymm_sbekk(arma::mat& r, const arma::mat& theta, int& max_iter, double& crit, arma::mat& signs) {

  arma::vec steps = {9.9,9,8,7,6,5,4,3,2,1,0.5,0.25,0.1,0.01,0.005,
                     0.001,0.0005,0.0001,0.00005,0.00001,0};
  double step = 0.1;
  int count_loop = 0;
  arma::mat theta_candidate = theta;
  int exit_loop = 0;
  arma::vec lik_all(max_iter+1, arma::fill::zeros);
  lik_all(0) = loglike_asymm_sbekk(theta, r, signs);


  while (count_loop < max_iter && exit_loop == 0) {
    arma::mat theta_loop = theta_candidate;
    arma::mat theta_temp = arma::zeros(theta_loop.n_rows, steps.n_elem);

    arma::mat score_function = score_asymm_sbekk(theta_loop, r, signs);
    arma::mat outer_score = score_function.t() * score_function;

    arma::mat outer_score_inv = arma::inv(outer_score);
    arma::mat score_function_sum = arma::sum(score_function);

    double lik = loglike_asymm_sbekk(theta_loop, r, signs);

    for (int i = 0; i < steps.n_elem; i++) {
      arma::vec temp = theta_candidate + step * steps(i) * outer_score_inv * score_function_sum.t();
      theta_temp.col(i) = temp;
    }


    arma::vec likelihood_candidates(steps.n_elem, arma::fill::zeros);
    likelihood_candidates(steps.n_elem - 1) = lik;

    int  j = steps.n_elem - 2;
    int exit_inner = 0;
    while (j >= 0 && exit_inner == 0) {
      likelihood_candidates(j) = loglike_asymm_sbekk(theta_temp.col(j), r, signs);
      //if (likelihood_candidates(j+1) > likelihood_candidates(j)) {
      //  exit_inner = 1;
      //}
      j -= 1;
    }

    //return likelihood_candidates;

    //int max_index = arma::index_max(likelihood_candidates.subvec(j, (steps.n_elem -1))) + j;
    int max_index = arma::index_max(likelihood_candidates);
    //return max_index;
    double likelihood_best = likelihood_candidates(max_index);

    // exit criterion strange
    //if (pow(likelihood_best - likelihood_candidates(steps.n_elem -1), 2)/abs(likelihood_candidates(steps.n_elem -1)) < crit) {
    //  exit_loop = 1;
    //}
    if (likelihood_best < lik_all(count_loop)) {
      exit_loop = 1;
      count_loop += 1;
    } else if (pow(likelihood_best - likelihood_candidates(steps.n_elem -1), 2)/abs(likelihood_candidates(steps.n_elem -1)) < crit) {// if (pow(likelihood_best - lik_all(count_loop), 2)/abs(lik_all(count_loop)) < crit) {
      exit_loop = 1;
      count_loop += 1;
      theta_candidate = theta_temp.col(max_index);
      lik_all(count_loop) = likelihood_candidates(steps.n_elem -1);
    } else {
      theta_candidate = theta_temp.col(max_index);
      count_loop += 1;
      lik_all(count_loop) = likelihood_candidates(steps.n_elem -1);
    }
  }

  double likelihood_final = loglike_asymm_sbekk(theta_candidate, r, signs);
  arma::mat score_final = score_asymm_sbekk(theta_candidate, r, signs);
  arma::mat s1_temp = arma::diagmat(arma::inv(score_final.t() * score_final));
  arma::mat s1 = arma::sqrt(s1_temp.diag());

  arma::mat t_val = theta_candidate/s1;
  return Rcpp::List::create(Rcpp::Named("theta") = theta_candidate,
                            Rcpp::Named("t_val") = t_val,
                            Rcpp::Named("likelihood") = likelihood_final,
                            Rcpp::Named("iter") = count_loop,
                            Rcpp::Named("likelihood_iter") = lik_all);
}

// [[Rcpp::export]]
arma::mat hesse_sbekk(arma::mat theta, arma::mat r){
  int n = r.n_rows;
  int N = r.n_cols;
  int N2 = pow(N,2);
  int NoOfVars_C = N*(N+1)/2;

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);
  arma::mat K_commutation = commutation_mat(N);
  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();


  arma::mat C1=arma::kron(arma::kron(arma::eye(N,N),K_commutation),arma::eye(N,N));
  arma::mat C2=2*arma::kron(arma::eye(N2,N2),D_duplication*D_gen_inv);
  arma::mat C3=C2*C1;

  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model
  arma::mat C = arma::zeros(N, N);
  int index =0;
  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }

  double a = theta[NoOfVars_C];
  double g = theta[NoOfVars_C+1];
  arma::mat CC = C * C.t();


  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  arma::mat dHda = arma::zeros(N2, 1);
  arma::mat dHdg = arma::zeros(N2, 1);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N,N)) * L_elimination.t();

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

  arma::mat ht_inv = arma::inv(ht);
  //Hessian
  arma::mat hessian = arma::zeros(theta.n_rows,theta.n_rows);

  //Second derivatives for t=1
  arma::mat dHdada = arma::zeros(N2, 1);
  arma::mat dHdadc = arma::zeros(N2, NoOfVars_C);
  arma::mat dHdadg = arma::zeros(N2, 1);

  arma::mat dHdgdg = arma::zeros(N2, 1);
  arma::mat dHdgdc = arma::zeros(N2, NoOfVars_C);
  arma::mat dHdgda = arma::zeros(N2, 1);

  arma::mat dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t();

  arma::mat dHdcda = arma::zeros(NoOfVars_C*N2,1);
  arma::mat dHdcdg = arma::zeros(NoOfVars_C*N2,1);


  arma::mat dHHdc=arma::join_horiz(dHdcdc,dHdcda,dHdcdg);

  arma::mat dHHda=arma::join_horiz(dHdadc,dHdada,dHdadg);

  arma::mat dHHdg=arma::join_horiz(dHdgdc,dHdgda,dHdgdg);


  arma::mat dHHdtheta=arma::join_vert(dHHdc,dHHda,dHHdg);

  arma::mat matt= arma::zeros(theta.n_rows,N2);
  matt(0,0)=1;
  arma::mat mat;

  for(int i=1; i < theta.n_rows; i++){
    mat=arma::zeros(theta.n_rows,N2);
    mat(i,0)=1;
    matt=arma::join_horiz(matt,mat);
  }
  arma::mat dHH=matt*dHHdtheta;

  for (int i=1; i<N2; i++){
    arma::mat matt2 = arma::join_horiz(arma::zeros(mat.n_rows, 1), matt.cols(0,matt.n_cols-2));
    dHH=arma::join_vert(dHH,matt2*dHHdtheta);
  }

  for(int i=0; i<theta.n_rows;i++){
    for(int j=0; j<theta.n_rows;j++){
      arma::mat dhi = dHdtheta.row(i).t();
      dhi = arma::reshape(dhi,N,N);

      arma::mat dhj = dHdtheta.row(j).t();
      dhj = arma::reshape(dhj,N,N);

      arma::mat temp = arma::zeros(1,N2);
      temp(0,0)=dHH(i,j);

      for(int k=1; k<N2; k++) {
        temp(0,k) = dHH(k*(theta.n_rows)+i,j);
      }

      temp= arma::reshape(temp,N,N);
      arma::mat ddh = temp.t();

      arma::mat mat1 = ddh*ht_inv-dhi*ht_inv*dhj*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dhj*ht_inv*dhi*ht_inv-r.row(0).t()*r.row(0)*ht_inv*ddh*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dhi*ht_inv*dhj*ht_inv;
      //update hessian
      hessian(i,j) = (-0.5)*arma::sum(mat1.diag());
    }
  }
  //Second derivatives for t>1

  for (int i=1; i<n;i++){

    dHdada = g * dHdada;
    dHdadg = g * dHdadg + dHda;
      //dHdadc = dHdadc; // Always zero (row may be deleted later on)

    dHdgda = dHdadg;
    dHdgdg = 2 * dHdg + g * dHdgdg;

    //dHdcda = dHdcda;// Always zero (row may be deleted later on)
    dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t()+g*dHdcdc;

    dHdcdg = arma::vectorise(dHdc)  + g*dHdcdg ;

    dHdgdc = arma::reshape(dHdcdg,N2, NoOfVars_C);

    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t() + g * dHdc;

    dHda = arma::reshape( r.row(i-1).t() * r.row(i-1),N2, 1) + g * dHda;

    dHdg = arma::reshape(ht,N2, 1) + g * dHdg;


    //update h
    ht = CC + a * r.row(i-1).t()*r.row(i-1)  + g * ht;

    //Computation of Hessian for each t>1


    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();


    arma::mat ht_inv=arma::inv(ht);

    arma::mat dHHdc=arma::join_horiz(dHdcdc,dHdcda,dHdcdg);
    arma::mat dHHda=arma::join_horiz(dHdadc,dHdada,dHdadg);
    arma::mat dHHdg=arma::join_horiz(dHdgdc,dHdgda,dHdgdg);

    arma::mat dHHdtheta=arma::join_vert(dHHdc,dHHda,dHHdg);


    arma::mat matt=arma::zeros(theta.n_rows,N2);
    matt(0,0)=1;
    arma::mat mat1;

    for (int j=1;j<theta.n_rows;j++){
      mat1=arma::zeros(theta.n_rows,N2);
      mat1(j,0)=1;
      matt=arma::join_horiz(matt,mat1);
    }
    arma::mat dHH=matt*dHHdtheta;

    for (int j=1;j<N2;j++){
      matt=arma::join_horiz(arma::zeros(mat1.n_rows,1),matt.cols(0,matt.n_cols-2));
      dHH=arma::join_vert(dHH,matt*dHHdtheta);
    }



    for (int l=0; l<theta.n_rows;l++){
      for (int j=0; j<theta.n_rows;j++){
        arma::mat dhi = dHdtheta.row(l).t();
        dhi = arma::reshape(dhi,N,N);

        arma::mat dhj = dHdtheta.row(j).t();
        dhj = arma::reshape(dhj,N,N);

        arma::mat temp = arma::zeros(1,N2);
        temp(0,0) = dHH(l,j);

        for (int k=1; k<N2; k++) {
          temp(0,k) = dHH(k*(theta.n_rows)+l,j);
        }

        temp = arma::reshape(temp,N,N);

        arma::mat ddh = temp.t();
        //get partial derivtatives for ll
        arma::mat mat = ddh*ht_inv-dhi*ht_inv*dhj*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dhj*ht_inv*dhi*ht_inv- r.row(i).t()*r.row(i)*ht_inv*ddh*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dhi*ht_inv*dhj*ht_inv;
        //update hessian
        hessian(l,j) = hessian(l,j)-0.5*arma::sum(mat.diag());

      }
    }

  }

  return hessian*(-1);

}

// [[Rcpp::export]]
arma::mat hesse_asymm_sbekk(arma::mat theta, arma::mat r, arma::mat signs){
  //int N = r.n_cols;
  //int N2 = pow(N,2);
  //int NoOfVars_C =2;
  double expected_indicator = expected_indicator_value(r, signs);

  // Partial derivatives for initial period t = 1
  int n = r.n_rows;
  int N = r.n_cols;
  int N2 = pow(N,2);
  int NoOfVars_C = N*(N+1)/2;

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);
  arma::mat K_commutation = commutation_mat(N);
  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();


  arma::mat C1=arma::kron(arma::kron(arma::eye(N,N),K_commutation),arma::eye(N,N));
  arma::mat C2=2*arma::kron(arma::eye(N2,N2),D_duplication*D_gen_inv);
  arma::mat C3=C2*C1;

  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model
  arma::mat C = arma::zeros(N, N);
  int index =0;
  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }

  double a = theta[NoOfVars_C];
  double b = theta[NoOfVars_C+1];
  double g = theta[NoOfVars_C+2];
  arma::mat CC = C * C.t();


  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  arma::mat dHda = arma::zeros(N2, 1);
  arma::mat dHdb = arma::zeros(N2, 1);
  arma::mat dHdg = arma::zeros(N2, 1);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N,N)) * L_elimination.t();

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdc, dHdg).t();

  arma::mat ht_inv = arma::inv(ht);
  //Hessian
  arma::mat hessian = arma::zeros(theta.n_rows,theta.n_rows);

  //Second derivatives for t=1
  arma::mat dHdada = arma::zeros(N2, 1);
  arma::mat dHdadb = arma::zeros(N2, 1);
  arma::mat dHdadc = arma::zeros(N2, NoOfVars_C);
  arma::mat dHdadg = arma::zeros(N2, 1);

  arma::mat dHdbda = arma::zeros(N2, 1);
  arma::mat dHdbdb = arma::zeros(N2, 1);
  arma::mat dHdbdc = arma::zeros(N2, NoOfVars_C);
  arma::mat dHdbdg = arma::zeros(N2, 1);

  arma::mat dHdgdg = arma::zeros(N2, 1);
  arma::mat dHdgdc = arma::zeros(N2, NoOfVars_C);
  arma::mat dHdgda = arma::zeros(N2, 1);
  arma::mat dHdgdb = arma::zeros(N2, 1);

  arma::mat dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t();
  arma::mat dHdcda = arma::zeros(NoOfVars_C*N2,1);
  arma::mat dHdcdb = arma::zeros(NoOfVars_C*N2,1);
  arma::mat dHdcdg = arma::zeros(NoOfVars_C*N2,1);


  arma::mat dHHdc=arma::join_horiz(dHdcdc,dHdcda,dHdcdb,dHdcdg);
  arma::mat dHHda=arma::join_horiz(dHdadc,dHdada,dHdadb,dHdadg);
  arma::mat dHHdb=arma::join_horiz(dHdbdc,dHdbda,dHdbdb,dHdbdg);
  arma::mat dHHdg=arma::join_horiz(dHdgdc,dHdgda,dHdgdb,dHdgdg);

  arma::mat dHHdtheta=arma::join_vert(dHHdc,dHHda,dHHdb,dHHdg);


  arma::mat matt= arma::zeros(theta.n_rows,N2);
  matt(0,0)=1;
  arma::mat mat;

  for(int i=1; i < theta.n_rows; i++){
    mat=arma::zeros(theta.n_rows,N2);
    mat(i,0)=1;
    matt=arma::join_horiz(matt,mat);
  }
  arma::mat dHH=matt*dHHdtheta;

  for (int i=1; i<N2; i++){
    arma::mat matt2 = arma::join_horiz(arma::zeros(mat.n_rows, 1), matt.cols(0,matt.n_cols-2));
    dHH=arma::join_vert(dHH,matt2*dHHdtheta);
  }

  for(int i=0; i<theta.n_rows;i++){
    for(int j=0; j<theta.n_rows;j++){
      arma::mat dhi = dHdtheta.row(i).t();
      dhi = arma::reshape(dhi,N,N);

      arma::mat dhj = dHdtheta.row(j).t();
      dhj = arma::reshape(dhj,N,N);

      arma::mat temp = arma::zeros(1,N2);
      temp(0,0)=dHH(i,j);

      for(int k=1; k<N2; k++) {
        temp(0,k) = dHH(k*(theta.n_rows)+i,j);
      }

      temp= arma::reshape(temp,N,N);
      arma::mat ddh = temp.t();

      arma::mat mat1 = ddh*ht_inv-dhi*ht_inv*dhj*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dhj*ht_inv*dhi*ht_inv-r.row(0).t()*r.row(0)*ht_inv*ddh*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dhi*ht_inv*dhj*ht_inv;
      //update hessian
      hessian(i,j) = (-0.5)*arma::sum(mat1.diag());
    }
  }
  //Second derivatives for t>1

  for (int i=1; i<n;i++){

    dHdbdg = g * dHdbdg + dHdb;
    dHdgdb = dHdbdg;

    dHdada = g * dHdada;
    dHdadg = g * dHdadg + dHda;
    //dHdadc = dHdadc; // Always zero (row may be deleted later on)

    dHdgda = dHdadg;
    dHdgdg = 2 * dHdg + g * dHdgdg;

    //dHdcda = dHdcda;// Always zero (row may be deleted later on)
    dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t()+g*dHdcdc;

    dHdcdg = arma::vectorise(dHdc)  + g*dHdcdg ;

    dHdgdc = arma::reshape(dHdcdg,N2, NoOfVars_C);

    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(C, arma::eye(N, N)) * L_elimination.t() + g * dHdc;

    dHda = arma::reshape( r.row(i-1).t() * r.row(i-1),N2, 1) + g * dHda;
    dHdb = indicatorFunction(r.row(i-1),signs) * arma::reshape( r.row(i-1).t() * r.row(i-1),N2, 1) + g * dHdb;

    dHdg = arma::reshape(ht,N2, 1) + g * dHdg;





    //dHdadc = dHdadc; // Always zero (row may be deleted later on)



    //dHdcda = dHdcda;// Always zero (row may be deleted later on)







    //update h
    ht = CC + (a+b*indicatorFunction(r.row(i-1),signs) ) * r.row(i-1).t()*r.row(i-1)  + g * ht;

    //Computation of Hessian for each t>1


    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdb, dHdg).t();


    arma::mat ht_inv=arma::inv(ht);

    arma::mat dHHdc=arma::join_horiz(dHdcdc,dHdcda,dHdcdb,dHdcdg);
    arma::mat dHHda=arma::join_horiz(dHdadc,dHdada,dHdadb,dHdadg);
    arma::mat dHHdb=arma::join_horiz(dHdbdc,dHdbda,dHdbdb,dHdbdg);
    arma::mat dHHdg=arma::join_horiz(dHdgdc,dHdgda,dHdgdb,dHdgdg);
    arma::mat dHHdtheta=arma::join_vert(dHHdc,dHHda,dHHdb,dHHdg);


    arma::mat matt=arma::zeros(theta.n_rows,N2);
    matt(0,0)=1;
    arma::mat mat1;

    for (int j=1;j<theta.n_rows;j++){
      mat1=arma::zeros(theta.n_rows,N2);
      mat1(j,0)=1;
      matt=arma::join_horiz(matt,mat1);
    }
    arma::mat dHH=matt*dHHdtheta;

    for (int j=1;j<N2;j++){
      matt=arma::join_horiz(arma::zeros(mat1.n_rows,1),matt.cols(0,matt.n_cols-2));
      dHH=arma::join_vert(dHH,matt*dHHdtheta);
    }



    for (int l=0; l<theta.n_rows;l++){
      for (int j=0; j<theta.n_rows;j++){
        arma::mat dhi = dHdtheta.row(l).t();
        dhi = arma::reshape(dhi,N,N);

        arma::mat dhj = dHdtheta.row(j).t();
        dhj = arma::reshape(dhj,N,N);

        arma::mat temp = arma::zeros(1,N2);
        temp(0,0) = dHH(l,j);

        for (int k=1; k<N2; k++) {
          temp(0,k) = dHH(k*(theta.n_rows)+l,j);
        }

        temp = arma::reshape(temp,N,N);

        arma::mat ddh = temp.t();
        //get partial derivtatives for ll
        arma::mat mat = ddh*ht_inv-dhi*ht_inv*dhj*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dhj*ht_inv*dhi*ht_inv- r.row(i).t()*r.row(i)*ht_inv*ddh*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dhi*ht_inv*dhj*ht_inv;
        //update hessian
        hessian(l,j) = hessian(l,j)-0.5*arma::sum(mat.diag());

      }
    }

  }

  return hessian*(-1);

}


// [[Rcpp::export]]
Rcpp::List sigma_sbekk(arma::mat& r, arma::mat& C, double a, double g) {
  // Computation of second order moment time paths and GARCH innovations
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat sigma = arma::zeros(r.n_rows, N2);
  arma::mat et = arma::zeros(r.n_rows, N);
  arma::mat omega = r.t() * r / r.n_rows;
  arma::mat ht = omega;
  sigma.row(0) = arma::vectorise(ht).t();

  //et.row(0) <- inv_gen(arma::sqrtmat_sympd(ht)) *  r.row(0).t();

  arma::mat CC  = C.t() * C;

  for (int i = 1; i < r.n_rows; i++) {
    ht = CC + a * r.row(i - 1).t() * r.row(i - 1)  + g*ht;
    sigma.row(i) = arma::vectorise(ht).t();
    et.row(i) = (inv_gen(arma::chol(ht).t()) *  r.row(i).t()).t();
  }

  return Rcpp::List::create(Rcpp::Named("sigma_t")= sigma,
                            Rcpp::Named("e_t") = et);
}

// [[Rcpp::export]]
Rcpp::List sigma_sbekk_asymm(arma::mat& r, arma::mat& C, double a, double b, double g, arma::mat signs) {
  // Computation of second order moment time paths and GARCH innovations
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat sigma = arma::zeros(r.n_rows, N2);
  arma::mat et = arma::zeros(r.n_rows, N);
  arma::mat omega = r.t() * r / r.n_rows;
  arma::mat ht = omega;
  sigma.row(0) = arma::vectorise(ht).t();
  double exp_signs = expected_indicator_value(r, signs);

  //et.row(0) <- inv_gen(arma::sqrtmat_sympd(ht)) * r.row(0).t();
  //changed CC to C.t() * C instead of C * C.t() because C is upper triagular in the inputs
  arma::mat CC  = C.t() * C;
  for (int i = 1; i < r.n_rows; i++) {
    ht = CC + (a+indicatorFunction(r.row(i-1),signs)*b) * r.row(i - 1).t() * r.row(i - 1)+ g * ht;
    sigma.row(i) = arma::vectorise(ht).t();
    et.row(i) = (inv_gen(arma::chol(ht).t()) * r.row(i).t()).t();
  }

  return Rcpp::List::create(Rcpp::Named("sigma_t")= sigma,
                            Rcpp::Named("e_t") = et);
}

//[[Rcpp::export]]
Rcpp::List random_grid_search_sBEKK(arma::mat r) {
  int n =r.n_cols;
  int N =r.n_rows;
  int l=0;
  int m=0;


  arma::mat C = arma::zeros(n,n);
  double a;
  double g;

  int numb_of_vars=2+n*(n+1)/2;

  arma::vec theta = arma::zeros(numb_of_vars,1);
  arma::vec thetaOptim=theta;
  arma::vec theta_mu = theta;
  //arma::vec theta_mu=theta;
  int counter= 0;
  int diagonal_elements = n;
  int diagonal_counter = 0;
  arma::mat uncond_var = r* r.t()/N;
  //set the initial expected values of the parameters
  for (int j=0; j < (n*(n+1)/2);j++){

    if(j == counter){
      theta_mu[j]=  0.05*uncond_var(j,j);
      counter+=diagonal_elements;
      diagonal_elements--;
    }

  }
  diagonal_counter=0;

  theta_mu[numb_of_vars-2]=0.2;
  theta_mu[numb_of_vars-1]=0.7;



  double best_val = loglike_sbekk(theta_mu,r);
  thetaOptim=theta_mu;
  //set the seed

  // Generating random values for A, B, C and G
  while(l<3000 && m<=17){
    int counter= 0;
    int diagonal_elements = n;
    int diagonal_counter = 0;

    for (int j=0; j < (n*(n+1)/2);j++){

      if(j == counter){
        theta[j]=  theta_mu[j]+arma::randn()*0.001;
        counter+=diagonal_elements;
        diagonal_elements--;
      }
      else{
        theta[j]=arma::randn()*0.00001+theta_mu[j];

      }
    }

    theta[numb_of_vars-2]= arma::randn()*0.03+theta_mu[numb_of_vars-2];

    theta[numb_of_vars-1]=arma::randn()*0.03+theta_mu[numb_of_vars-1];




    //arma::mat C = arma::zeros(n,n);
    int  index=0;
    for(int j=0; j <n; j++){
      for (int k = j;k < n; k++) {
        C(k,j)=theta[index];
        index++;
      }
    }

    a = theta[numb_of_vars-2];
    g = theta[numb_of_vars-1];

    if(valid_sbekk(C,a,g)){
      l++;
      double llv=loglike_sbekk(theta,r);

      if(llv>best_val){

        m++;
        best_val=llv;
        thetaOptim=theta;
        theta_mu=thetaOptim;

      }
      if(l>500 || m>=5){
        theta_mu=thetaOptim;
      }
    }

  }
  return Rcpp::List::create(Rcpp::Named("thetaOptim") = thetaOptim,
                            Rcpp::Named("best_val") = best_val);

}

//[[Rcpp::export]]
Rcpp::List random_grid_search_asymmetric_sBEKK(arma::mat r, arma::mat signs) {
  int n =r.n_cols;
  int N =r.n_rows;
  int l=0;
  int m=0;


  arma::mat C = arma::zeros(n,n);
  double a;
  double b;
  double g;

  int numb_of_vars=3+n*(n+1)/2;

  arma::vec theta = arma::zeros(numb_of_vars,1);
  arma::vec thetaOptim=theta;
  arma::vec theta_mu = theta;
  //arma::vec theta_mu=theta;
  int counter= 0;
  int diagonal_elements = n;
  int diagonal_counter = 0;
  arma::mat uncond_var = r* r.t()/N;
  //set the initial expected values of the parameters
  for (int j=0; j < (n*(n+1)/2);j++){

    if(j == counter){
      theta_mu[j]=  0.05*uncond_var(j,j);
      counter+=diagonal_elements;
      diagonal_elements--;

    }
  }
      theta_mu[numb_of_vars-3]=0.1;


      theta_mu[numb_of_vars-2]=0.1;

      theta_mu[numb_of_vars-1]=0.7;



      double best_val = loglike_asymm_sbekk(theta_mu,r,signs);
      thetaOptim=theta_mu;
      //set the seed

      // Generating random values for A, B, C and G
      while(l<4000 && m<=17){
        int counter= 0;
        int diagonal_elements = n;
        int diagonal_counter = 0;

        for (int j=0; j < (n*(n+1)/2);j++){

          if(j == counter){
            theta[j]=  theta_mu[j]+arma::randn()*0.001;
            counter+=diagonal_elements;
            diagonal_elements--;
          }
          else{
            theta[j]=arma::randn()*0.00001+theta_mu[j];

          }
        }

        theta[numb_of_vars-3]= arma::randn()*0.03+theta_mu[numb_of_vars-3];

        theta[numb_of_vars-2]= arma::randn()*0.03+theta_mu[numb_of_vars-2];

        theta[numb_of_vars-1]=arma::randn()*0.03+theta_mu[numb_of_vars-1];




        //arma::mat C = arma::zeros(n,n);
        int  index=0;
        for(int j=0; j <n; j++){
          for (int k = j;k < n; k++) {
            C(k,j)=theta[index];
            index++;
          }
        }

        a = theta[numb_of_vars-3];
        b = theta[numb_of_vars-2];
        g = theta[numb_of_vars-1];

        if(valid_asymm_sbekk(C,a,b,g,r,signs)){
          l++;
          double llv=loglike_asymm_sbekk(theta,r,signs);

          if(llv>best_val){
            m++;
            best_val=llv;
            thetaOptim=theta;
            theta_mu=thetaOptim;
          }
          if(l>100 || m>=5){
            theta_mu=thetaOptim;
          }
        }
      }
        return Rcpp::List::create(Rcpp::Named("thetaOptim") = thetaOptim,
                                  Rcpp::Named("best_val") = best_val);

      }
