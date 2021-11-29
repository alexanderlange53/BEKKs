#include <RcppArmadillo.h>
#include "IndicatorFunctions.h"




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
bool valid_scalar_bekk(double a, double b){
  double sum=a+b;
  if(sum >= 1){
    return false;
  }


  if(a<=0 || b<=0) {
    return false;
  }
  else{
    return true;
  }
}

// [[Rcpp::export]]
double loglike_scalar_bekk(const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2;
  double a = theta[0];
  double b = theta[1];


  // check constraints
  if (valid_scalar_bekk(a,b) == false) {
    return -1e25;
  }

  // compute inital H
  arma::mat  sample_covariance = (r.t() * r) / r.n_rows;
  arma::mat  H = sample_covariance;


  double llv = arma::as_scalar(log(arma::det(H)) + r.row(0) * arma::inv(H) * r.row(0).t());
  for (int i = 1; i < NoOBs; i++) {
    H = (1-a-b)*sample_covariance + a * r.row(i - 1).t() * r.row(i - 1)  + b * H ;
    llv += arma::as_scalar(log(arma::det(H)) + r.row(i) * arma::inv(H) * r.row(i).t());
  }


  return -0.5 * n * NoOBs * log(2 * M_PI) - 0.5 * llv;
}

// [[Rcpp::export]]
arma::mat score_scalar_bekk(const arma::mat& theta, arma::mat& r) {

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);


  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model

  // Partial derivatives for initial period t = 1
  arma::mat sample_covar = r.t() * r / r.n_rows;
  arma::mat ht = sample_covar;

  arma::mat dHda = 0;

  arma::mat dHdb = 0;

  double a = theta[0];
  double b = theta[1];



    arma::mat ht_sqrt_inv = arma::inv(ht);
  //arma::vec et = ht_sqrt_inv * r.row(0).t();

  //return dHdtheta;
  // Score function

    //Rcpp::Rcout << ht_sqrt_inv;
    //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
    arma::mat mat_temp = dHda * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dHda * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, 0) = -(0.5) * arma::sum(mat_temp.diag());
    mat_temp = dHdb * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dHdb * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, 1) = -(0.5) * arma::sum(mat_temp.diag());


  // Partial derivatives for period t >= 2

  for (int i = 1; i < r.n_rows; i++) {
    dHda = r.row(i-1).t() * r.row(i-1)-sample_covar+b*dHda;



    dHdb = ht-sample_covar+b*dHdb;

    //dHda =  arma::kron(arma::eye(N, N),at*r.row(i-1).t() *r.row(i-1)) + arma::kron(at * r.row(i-1).t() *r.row(i-1), arma::eye(N, N)) * commutation_mat(N)+ GGt * dHda;


    ht = (1-a-b)*sample_covar+a*r.row(i-1).t() * r.row(i-1)+b*ht;

    //ht_sqrt_inv = arma::inv(arma::real(arma::sqrtmat(ht)));
    arma::mat ht_sqrt_inv = arma::inv(ht);
    //et = ht_sqrt_inv * r.row(i).t();





      //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
      arma::mat mat_temp = dHda * ht_sqrt_inv - r.row(i).t() * r.row(i)* ht_sqrt_inv * dHda * ht_sqrt_inv;
      gradients(i, 0) = -(0.5) * arma::sum(mat_temp.diag());

      mat_temp = dHdb * ht_sqrt_inv - r.row(i).t() * r.row(i)* ht_sqrt_inv * dHdb * ht_sqrt_inv;
      gradients(i, 1) = -(0.5) * arma::sum(mat_temp.diag());
    }


  return gradients;
}

// [[Rcpp::export]]
Rcpp::List  bhh_scalar_bekk(arma::mat& r, const arma::mat& theta, int& max_iter, double& crit) {

  arma::vec steps = {9.9,9,8,7,6,5,4,3,2,1,0.5,0.25,0.1,0.01,0.005,
                     0.001,0.0005,0.0001,0.00005,0.00001,0};
  double step = 0.1;
  int count_loop = 0;
  arma::mat theta_candidate = theta;
  int exit_loop = 0;
  arma::vec lik_all(max_iter+1, arma::fill::zeros);
  lik_all(0) = loglike_scalar_bekk(theta, r);


  while (count_loop < max_iter && exit_loop == 0) {
    arma::mat theta_loop = theta_candidate;
    arma::mat theta_temp = arma::zeros(theta_loop.n_rows, steps.n_elem);

    arma::mat score_function = score_scalar_bekk(theta_loop, r);
    arma::mat outer_score = score_function.t() * score_function;

    arma::mat outer_score_inv = arma::inv(outer_score);
    arma::mat score_function_sum = arma::sum(score_function);

    double lik = loglike_scalar_bekk(theta_loop, r);

    for (int i = 0; i < steps.n_elem; i++) {
      arma::vec temp = theta_candidate + step * steps(i) * outer_score_inv * score_function_sum.t();
      theta_temp.col(i) = temp;
    }


    arma::vec likelihood_candidates(steps.n_elem, arma::fill::zeros);
    likelihood_candidates(steps.n_elem - 1) = lik;

    int  j = steps.n_elem - 2;
    int exit_inner = 0;
    while (j >= 0 && exit_inner == 0) {
      likelihood_candidates(j) = loglike_scalar_bekk(theta_temp.col(j), r);
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

  double likelihood_final = loglike_scalar_bekk(theta_candidate, r);
  arma::mat score_final = score_scalar_bekk(theta_candidate, r);
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
arma::mat hesse_scalar_bekk(arma::mat theta, arma::mat r){
  int n = r.n_rows;
  int N = r.n_cols;
  int N2 = pow(N,2);
  int NoOfVars_C =2;


  // Partial derivatives for initial period t = 1
  arma::mat sample_covar = r.t() * r / r.n_rows;
  arma::mat ht = sample_covar;

  arma::mat dHda = 0;
  arma::mat dHdb = 0;

  double a = theta[0];
  double b = theta[1];


  arma::mat ht_inv = arma::inv(ht);
  //Hessian
  arma::mat hessian = arma::zeros(theta.n_rows,theta.n_rows);

  //Second derivatives for t=1
  arma::mat dHdada = 0;
  arma::mat dHdadb = 0;
  arma::mat dHdbdb = 0;



      arma::mat mat1 = dHdada*ht_inv-dHda*ht_inv*dHda*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHda*ht_inv*dHda*ht_inv-r.row(0).t()*r.row(0)*ht_inv*dHdada*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHda*ht_inv*dHda*ht_inv;
      //update hessian
      hessian(0,0) = (-0.5)*arma::sum(mat1.diag());

      mat1 = dHdbdb*ht_inv-dHdb*ht_inv*dHdb*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHdb*ht_inv*dHdb*ht_inv-r.row(0).t()*r.row(0)*ht_inv*dHdbdb*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHdb*ht_inv*dHdb*ht_inv;
      //update hessian
      hessian(1,1) = (-0.5)*arma::sum(mat1.diag());

      mat1 = dHdada*ht_inv-dHda*ht_inv*dHdb*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHdb*ht_inv*dHda*ht_inv-r.row(0).t()*r.row(0)*ht_inv*dHdada*ht_inv+r.row(0).t()*r.row(0)*ht_inv*dHda*ht_inv*dHdb*ht_inv;
      //update hessian
      hessian(0,1) = (-0.5)*arma::sum(mat1.diag());

      hessian(0,1)=hessian(1,0);


  //Second derivatives for t>1

  for (int i=1; i<n;i++){

    dHdada = b * dHdada;
    dHdbdb = 2*dHdb+b*dHdbdb;
    dHdadb = b * dHdadb +dHda ;
    //update h
    ht = (1-a-b)*sample_covar+ a * r.row(i-1).t()*r.row(i-1) + b * ht;

    //Computation of Hessian for each t>1





    arma::mat mat1 = dHdada*ht_inv-dHda*ht_inv*dHda*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHda*ht_inv*dHda*ht_inv-r.row(i).t()*r.row(i)*ht_inv*dHdada*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHda*ht_inv*dHda*ht_inv;
    //update hessian
    hessian(0,0) = hessian(0,0)-0.5*arma::sum(mat1.diag());

    mat1 = dHdbdb*ht_inv-dHdb*ht_inv*dHdb*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHdb*ht_inv*dHdb*ht_inv-r.row(i).t()*r.row(i)*ht_inv*dHdbdb*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHdb*ht_inv*dHdb*ht_inv;
    //update hessian
    hessian(1,1) = hessian(1,1)-0.5*arma::sum(mat1.diag());

    mat1 = dHdada*ht_inv-dHda*ht_inv*dHdb*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHdb*ht_inv*dHda*ht_inv-r.row(i).t()*r.row(i)*ht_inv*dHdada*ht_inv+r.row(i).t()*r.row(i)*ht_inv*dHda*ht_inv*dHdb*ht_inv;
    //update hessian
    hessian(0,1) = hessian(0,1) -0.5*arma::sum(mat1.diag());

    hessian(0,1)=hessian(1,0);  //update hessian





  }

  return hessian*(-1);

}
