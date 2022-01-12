#include <RcppArmadillo.h>

#include "IndicatorFunctions.h"



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
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

// [[Rcpp::export]]
bool valid_bekk(arma::mat& C,arma::mat& A,arma::mat& G){
  int n =C.n_cols;
  arma::mat prod = kron(A,A)+kron(G,G);

  arma::vec eigvals;
  eigvals= abs(arma::eig_gen(prod));
  double max=0;
  for (int i=0; i< eigvals.n_elem; i++){
    if(eigvals[i]>max){
      max=eigvals[i];
    }
  }
  if(max >= 1){
    return false;
  }

  for (int i=0; i<n;i++){
    if(C(i,i)<=0){
      return false;
    }
  }
  if(A(0,0)<=0 || G(0,0)<=0) {
    return false;
  }
  else{
    return true;
  }
}

// [[Rcpp::export]]
double expected_indicator_value(arma::mat r, arma::mat signs){
  int N =r.n_rows;
  double exp_indicator_value = 0;
  for (int i=0; i<N;i++){
    exp_indicator_value+=indicatorFunction(r.row(i),signs);
  }
  exp_indicator_value=exp_indicator_value/N;
  return exp_indicator_value;
}

// [[Rcpp::export]]
bool valid_asymm_bekk(arma::mat& C,arma::mat& A, arma::mat& B ,arma::mat& G, arma::mat r, arma::mat signs){
  int n =C.n_cols;
  int N =r.n_rows;
  double exp_indicator_value = expected_indicator_value(r,signs);

  arma::mat prod = kron(A,A)+exp_indicator_value*kron(B,B)+kron(G,G);

  arma::vec eigvals;
  eigvals= abs(arma::eig_gen(prod));
  double max=0;
  for (int i=0; i< eigvals.n_elem; i++){
    if(eigvals[i]>max){
      max=eigvals[i];
    }
  }

  if(max >= 1){
    return false;
  }

  for (int i=0; i<n;i++){
    if(C(i,i)<=0){
      return false;
    }
  }
  if(A(0,0)<=0 || B(0,0)<=0 || G(0,0)<=0) {
    return false;
  }

  else{
    return true;
  }
}


// [[Rcpp::export]]
bool valid_asymm_bekk_sim(arma::mat& C,arma::mat& A, arma::mat& B ,arma::mat& G, double exp_indicator_value, arma::mat signs){



  arma::mat prod = kron(A,A)+exp_indicator_value*kron(B,B)+kron(G,G);

  arma::vec eigvals;
  eigvals= abs(arma::eig_gen(prod));
  double max=0;
  for (int i=0; i< eigvals.n_elem; i++){
    if(eigvals[i]>max){
      max=eigvals[i];
    }
  }

  if(max >= 1){
    return false;
  }

  for (int i=0; i < C.n_cols;i++){
    if(C(i,i)<=0){
      return false;
    }
  }
  if(A(0,0)<=0 || B(0,0)<=0 || G(0,0)<=0) {
    return false;
  }

  else{
    return true;
  }
}

// [[Rcpp::export]]
double loglike_bekk(const arma::vec& theta, const arma::mat& r) {
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


   arma::mat A = arma::reshape(theta.subvec(index, (index + pow(n, 2)) - 1 ).t(), n, n);
   arma::mat G = arma::reshape(theta.subvec((index +  pow(n, 2)), numb_of_vars-1).t(), n, n);

// check constraints
    if (valid_bekk(C, A, G) == false) {
      return -1e25;
    }

// compute inital H
    arma::mat H = (r.t() * r) / r.n_rows;

    arma::mat CC  = C * C.t();
    arma::mat At  = A.t();
    arma::mat Gt  = G.t();

    double llv = arma::as_scalar(log(arma::det(H)) + r.row(0) * inv_gen(H) * r.row(0).t());
    for (int i = 1; i < NoOBs; i++) {
      H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
      llv += arma::as_scalar(log(arma::det(H)) + r.row(i) * inv_gen(H) * r.row(i).t());
    }


    return -0.5 * n * NoOBs * log(2 * M_PI) - 0.5 * llv;
}

// [[Rcpp::export]]
double loglike_asymm_bekk(const arma::vec& theta, const arma::mat& r, arma::mat& signs) {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 3 * pow(n, 2) + n * (n + 1)/2;
  arma::mat C = arma::zeros(n, n);

  int index = 0;

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      C(j, i) = theta(index);
      index += 1;
    }
  }


  arma::mat A = arma::reshape(theta.subvec(index, (index + pow(n, 2)) - 1 ).t(), n, n);
  arma::mat B = arma::reshape(theta.subvec(index +  pow(n, 2), index +  2*pow(n, 2) -1).t(), n, n);
  arma::mat G = arma::reshape(theta.subvec(index +  2*pow(n, 2), numb_of_vars-1).t(), n, n);

  // check constraints
  if (valid_asymm_bekk(C, A, B , G, r, signs) == false) {
    return -1e25;
  }
  //transform r to r star for asymmetry handling

  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;

  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Bt  = B.t();
  arma::mat Gt  = G.t();

  double llv = arma::as_scalar(log(arma::det(H)) + r.row(0) * inv_gen(H) * r.row(0).t());
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i - 1), signs) * Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * H * G;
    llv += arma::as_scalar(log(arma::det(H)) + r.row(i) * inv_gen(H) * r.row(i).t());
  }


  return -0.5 * n * NoOBs * log(2 * M_PI) - 0.5 * llv;
}

// [[Rcpp::export]]
arma::mat score_bekk(const arma::mat& theta, arma::mat& r) {
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);


  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);


  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model
  arma::mat c0 = arma::zeros(N, N);
  //arma::uvec lw_idx = arma::trimatu_ind(arma::size(c0));
  //c0.elem(lw_idx) = theta.rows(0, (N * (N + 1)/2) - 1);
  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      c0(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat a = arma::reshape(theta.rows((N * (N+1)/2), (pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat g = arma::reshape(theta.rows(((pow(N, 2) + (N * (N + 1)/2))), (2*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  c0=c0.t();
  arma::mat c_full = c0.t() * c0;
  arma::mat at = a.t();
  arma::mat gt = g.t();
  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  //arma::mat dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), a.t() * ht);
  arma::mat dHda = arma::zeros(N2, N2);

  arma::mat dHdg = arma::zeros(N2, N2);
  // arma::mat dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), g.t() * ht);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N, N)) * L_elimination.t();


  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

  //arma::mat ht_sqrt_inv = arma::inv(arma::real(arma::sqrtmat(ht)));
  arma::mat ht_sqrt_inv = inv_gen(ht);
  //arma::vec et = ht_sqrt_inv * r.row(0).t();

  //return dHdtheta;
  // Score function
  for (int k = 0; k < theta.n_rows; k++) {

    arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
    //Rcpp::Rcout << ht_sqrt_inv;
    //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
    arma::mat mat_temp = dh * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dh * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, k) = -(0.5) * arma::sum(mat_temp.diag());
    //Rcpp::Rcout << gradients(0, k);
  }
  //Rcpp::Rcout << "Matrix M\n" << gradients;
  // Partial derivatives for period t >= 2
  arma::mat GGt = arma::kron(g, g).t();
    for (int i = 1; i < r.n_rows; i++) {
      dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N, N)) * L_elimination.t() + GGt * dHdc;



    dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), a.t() * r.row(i-1).t() * r.row(i-1)) + GGt * dHda;

    //dHda =  arma::kron(arma::eye(N, N),at*r.row(i-1).t() *r.row(i-1)) + arma::kron(at * r.row(i-1).t() *r.row(i-1), arma::eye(N, N)) * commutation_mat(N)+ GGt * dHda;

   dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), g.t() * ht) + GGt * dHdg;

    //dHdg = arma::kron(arma::eye(N, N),gt*ht) + arma::kron(gt * ht , arma::eye(N, N)) * commutation_mat(N) + GGt * dHdg;

    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

    ht = c_full + a.t() * r.row(i-1).t() * r.row(i-1) * a + g.t() * ht * g;

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
arma::mat score_asymm_bekk(const arma::mat& theta, arma::mat& r, arma::mat& signs) {
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);


  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);


  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model
  arma::mat c0 = arma::zeros(N, N);
  //arma::uvec lw_idx = arma::trimatu_ind(arma::size(c0));
  //c0.elem(lw_idx) = theta.rows(0, (N * (N + 1)/2) - 1);
  int index = 0;

  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      c0(j, i) = theta(index);
      index += 1;
    }
  }
  arma::mat a = arma::reshape(theta.rows((N * (N+1)/2), (pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat b = arma::reshape(theta.rows(pow(N, 2) + (N * (N + 1)/2), (2*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat g = arma::reshape(theta.rows(2*pow(N, 2) + (N * (N + 1)/2), (3*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);

  c0=c0.t();
  arma::mat c_full = c0.t() * c0;
  arma::mat at = a.t();
  arma::mat bt = b.t();
  arma::mat gt = g.t();
  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  //arma::mat dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), a.t() * ht);
  arma::mat dHda = arma::zeros(N2, N2);
  arma::mat dHdb = arma::zeros(N2, N2);
  arma::mat dHdg = arma::zeros(N2, N2);
  // arma::mat dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), g.t() * ht);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N, N)) * L_elimination.t();


  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda,dHdb, dHdg).t();

  //arma::mat ht_sqrt_inv = arma::inv(arma::real(arma::sqrtmat(ht)));
  arma::mat ht_sqrt_inv = inv_gen(ht);
  //arma::vec et = ht_sqrt_inv * r.row(0).t();

  //return dHdtheta;
  // Score function
  for (int k = 0; k < theta.n_rows; k++) {

    arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);
    //Rcpp::Rcout << ht_sqrt_inv;
    //arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());
    arma::mat mat_temp = dh * ht_sqrt_inv - r.row(0).t() * r.row(0) * ht_sqrt_inv * dh * ht_sqrt_inv;
    //Rcpp::Rcout << mat_temp;
    gradients(0, k) = -(0.5) * arma::sum(mat_temp.diag());
    //Rcpp::Rcout << gradients(0, k);
  }
  //Rcpp::Rcout << "Matrix M\n" << gradients;
  // Partial derivatives for period t >= 2
  arma::mat GGt = arma::kron(g, g).t();
  for (int i = 1; i < r.n_rows; i++) {
    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N, N)) * L_elimination.t() + GGt * dHdc;



    dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), at * r.row(i-1).t() * r.row(i-1)) + GGt * dHda;
    dHdb = indicatorFunction(r.row(i-1),signs) *  2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), bt * r.row(i-1).t() * r.row(i-1)) + GGt * dHdb;

    //dHda =  arma::kron(arma::eye(N, N),at*r.row(i-1).t() *r.row(i-1)) + arma::kron(at * r.row(i-1).t() *r.row(i-1), arma::eye(N, N)) * commutation_mat(N)+ GGt * dHda;

    dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), gt * ht) + GGt * dHdg;

    //dHdg = arma::kron(arma::eye(N, N),gt*ht) + arma::kron(gt * ht , arma::eye(N, N)) * commutation_mat(N) + GGt * dHdg;

    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda,dHdb, dHdg).t();

    ht = c_full + a.t() * r.row(i-1).t() * r.row(i-1) * a +indicatorFunction(r.row(i-1),signs)*b.t() * r.row(i-1).t() * r.row(i-1) * b+ g.t() * ht * g;

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
Rcpp::List  bhh_bekk(arma::mat& r, const arma::mat& theta, int& max_iter, double& crit) {

  arma::vec steps = {9.9,9,8,7,6,5,4,3,2,1,0.5,0.25,0.1,0.01,0.005,
                     0.001,0.0005,0.0001,0.00005,0.00001,0};
  double step = 0.1;
  int count_loop = 0;
  arma::mat theta_candidate = theta;
  int exit_loop = 0;
  arma::vec lik_all(max_iter+1, arma::fill::zeros);
  lik_all(0) = loglike_bekk(theta, r);


  while (count_loop < max_iter && exit_loop == 0) {
    arma::mat theta_loop = theta_candidate;
    arma::mat theta_temp = arma::zeros(theta_loop.n_rows, steps.n_elem);

    arma::mat score_function = score_bekk(theta_loop, r);
    arma::mat outer_score = score_function.t() * score_function;

    arma::mat outer_score_inv = inv_gen(outer_score);
    arma::mat score_function_sum = arma::sum(score_function);

    double lik = loglike_bekk(theta_loop, r);

     for (int i = 0; i < steps.n_elem; i++) {
        arma::vec temp = theta_candidate + step * steps(i) * outer_score_inv * score_function_sum.t();
        theta_temp.col(i) = temp;
      }


      arma::vec likelihood_candidates(steps.n_elem, arma::fill::zeros);
      likelihood_candidates(steps.n_elem - 1) = lik;

      int  j = steps.n_elem - 2;
      int exit_inner = 0;
      while (j >= 0 && exit_inner == 0) {
        likelihood_candidates(j) = loglike_bekk(theta_temp.col(j), r);
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

  double likelihood_final = loglike_bekk(theta_candidate, r);
  arma::mat score_final = score_bekk(theta_candidate, r);
  arma::mat s1_temp = arma::diagmat(inv_gen(score_final.t() * score_final));
  arma::mat s1 = arma::sqrt(s1_temp.diag());

  arma::mat t_val = theta_candidate/s1;
  return Rcpp::List::create(Rcpp::Named("theta") = theta_candidate,
                       Rcpp::Named("t_val") = t_val,
                       Rcpp::Named("likelihood") = likelihood_final,
                       Rcpp::Named("iter") = count_loop,
                       Rcpp::Named("likelihood_iter") = lik_all);
}

// [[Rcpp::export]]
Rcpp::List  bhh_asymm_bekk(arma::mat& r, const arma::mat& theta, int& max_iter, double& crit, arma::mat& signs) {

  arma::vec steps = {9.9,9,8,7,6,5,4,3,2,1,0.5,0.25,0.1,0.01,0.005,
                     0.001,0.0005,0.0001,0.00005,0.00001,0};
  double step = 0.1;
  int count_loop = 0;
  arma::mat theta_candidate = theta;
  int exit_loop = 0;
  arma::vec lik_all(max_iter+1, arma::fill::zeros);
  lik_all(0) = loglike_asymm_bekk(theta, r, signs);


  while (count_loop < max_iter && exit_loop == 0) {
    arma::mat theta_loop = theta_candidate;
    arma::mat theta_temp = arma::zeros(theta_loop.n_rows, steps.n_elem);

    arma::mat score_function = score_asymm_bekk(theta_loop, r, signs);
    arma::mat outer_score = score_function.t() * score_function;

    arma::mat outer_score_inv = inv_gen(outer_score);
    arma::mat score_function_sum = arma::sum(score_function);

    double lik = loglike_asymm_bekk(theta_loop, r, signs);

    for (int i = 0; i < steps.n_elem; i++) {
      arma::vec temp = theta_candidate + step * steps(i) * outer_score_inv * score_function_sum.t();
      theta_temp.col(i) = temp;
    }


    arma::vec likelihood_candidates(steps.n_elem, arma::fill::zeros);
    likelihood_candidates(steps.n_elem - 1) = lik;

    int  j = steps.n_elem - 2;
    int exit_inner = 0;
    while (j >= 0 && exit_inner == 0) {
      likelihood_candidates(j) = loglike_asymm_bekk(theta_temp.col(j), r, signs);
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

  double likelihood_final = loglike_asymm_bekk(theta_candidate, r, signs);
  arma::mat score_final = score_asymm_bekk(theta_candidate, r, signs);
  arma::mat s1_temp = arma::diagmat(inv_gen(score_final.t() * score_final));
  arma::mat s1 = arma::sqrt(s1_temp.diag());

  arma::mat t_val = theta_candidate/s1;
  return Rcpp::List::create(Rcpp::Named("theta") = theta_candidate,
                            Rcpp::Named("t_val") = t_val,
                            Rcpp::Named("likelihood") = likelihood_final,
                            Rcpp::Named("iter") = count_loop,
                            Rcpp::Named("likelihood_iter") = lik_all);
}



//[[Rcpp::export]]
Rcpp::List random_grid_search_BEKK(arma::mat r) {
  int n =r.n_cols;
  int N =r.n_rows;
  int l=0;
  int m=0;


  arma::mat C = arma::zeros(n,n);
  arma::mat A = arma::zeros(n,n);
  arma::mat G = arma::zeros(n,n);

  int numb_of_vars=2*(pow(n,2))+n*(n+1)/2;

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

  for (int j=(n*(n+1)/2); j < (numb_of_vars-2*n*n);j++){
    if(j == (n*(n+1)/2+diagonal_counter*(n+1))){
      diagonal_counter++;
      theta_mu[j]=0.3;
    }
    else{
      theta_mu[j]=0.001;
    }
  }
  diagonal_counter=0;
  for (int j=(n*(n+1)/2); j < (numb_of_vars-2*n*n);j++){
    if(j == (n*(n+1)/2+diagonal_counter*(n+1))){
      diagonal_counter++;
      theta_mu[j+n*n]=0.3;
    }
    else{
      theta_mu[j+n*n]=0.001;
    }
  }



  double best_val = loglike_bekk(theta_mu,r);
  thetaOptim=theta_mu;
  //set the seed

  // Generating random values for A, B, C and G
  while(l<10000 && m<=17){
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
    diagonal_counter=0;
    for (int j=(n*(n+1)/2); j < numb_of_vars;j++){
      if(j == (n*(n+1)/2+diagonal_counter*(n+1)) && j< ((n*(n+1)/2) +n*n)){
        diagonal_counter++;
        theta[j]= arma::randn()*0.001+theta_mu[j];
      }
      else if(j == (n*(n+1)/2+n*n+(diagonal_counter-n)*(n+1)) && j>=(n*(n+1)/2+n*n)){
        diagonal_counter++;
        theta[j]= arma::randn()*0.001+theta_mu[j];
      }
      else{
        theta[j]=arma::randn()*0.00001+theta_mu[j];
      }
    }


    //arma::mat C = arma::zeros(n,n);
    int  index=0;
    for(int j=0; j <n; j++){
      for (int k = j;k < n; k++) {
        C(k,j)=theta[index];
        index++;
      }
    }

    A = arma::reshape(theta.rows((n * (n+1)/2), pow(n, 2) + (n * (n + 1)/2) - 1), n, n);
    G = arma::reshape(theta.rows(pow(n, 2) + (n * (n + 1)/2), 2*pow(n, 2) + (n * (n + 1)/2) - 1), n, n);

    if(valid_bekk(C,A,G)){
      l++;
      double llv=loglike_bekk(theta,r);

      if(llv>best_val){

        m++;
        best_val=llv;
        thetaOptim=theta;
        theta_mu=thetaOptim;

      }
      if(l>2000 || m>=5){
        theta_mu=thetaOptim;
      }
    }

  }
  return Rcpp::List::create(Rcpp::Named("thetaOptim") = thetaOptim,
                            Rcpp::Named("best_val") = best_val);

}


//[[Rcpp::export]]
Rcpp::List random_grid_search_asymmetric_BEKK(arma::mat r, int nc, arma::mat signs) {
  int n =r.n_cols;
  int N =r.n_rows;
  int l=0;
  int m=0;


  arma::mat C = arma::zeros(n,n);
  arma::mat A = arma::zeros(n,n);
  arma::mat B = arma::zeros(n,n);
  arma::mat G = arma::zeros(n,n);

  int numb_of_vars=3*(pow(n,2))+n*(n+1)/2;

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

  for (int j=(n*(n+1)/2); j < (numb_of_vars-2*n*n);j++){
    if(j == (n*(n+1)/2+diagonal_counter*(n+1))){
      diagonal_counter++;
      theta_mu[j]=0.3;
    }
    else{
      theta_mu[j]=0.001;
    }
  }
  diagonal_counter=0;
  for (int j=(n*(n+1)/2); j < (numb_of_vars-2*n*n);j++){
    if(j == (n*(n+1)/2+diagonal_counter*(n+1))){
      diagonal_counter++;
      theta_mu[j+n*n]=0.3;
    }
    else{
      theta_mu[j+n*n]=0.001;
    }
  }
  diagonal_counter=0;
  for (int j=(n*(n+1)/2); j < (numb_of_vars-2*n*n);j++){
    if(j == (n*(n+1)/2+diagonal_counter*(n+1))){
      diagonal_counter++;
      theta_mu[j+n*n*2]=0.92;
    }
    else{
      theta_mu[j+n*n*2]=0.001;
    }
  }


  double best_val = loglike_asymm_bekk(theta_mu,r,signs);
  thetaOptim=theta_mu;
  //set the seed

  // Generating random values for A, B, C and G
  while(l<10000/log(2+nc) && m<=17){
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
    diagonal_counter=0;
    for (int j=(n*(n+1)/2); j < numb_of_vars;j++){
      if(j == (n*(n+1)/2+diagonal_counter*(n+1)) && j< ((n*(n+1)/2) +n*n)){
        diagonal_counter++;
        theta[j]= arma::randn()*0.0001+theta_mu[j];
      }
      else if(j == (n*(n+1)/2+n*n+(diagonal_counter-n)*(n+1)) && j>=(n*(n+1)/2+n*n) && j<((n*(n+1)/2) +2*n*n)){
        diagonal_counter++;
        theta[j]= arma::randn()*0.01+theta_mu[j];
      }
      else if(j == (n*(n+1)/2+2*n*n+(diagonal_counter-2*n)*(n+1)) && j>=(n*(n+1)/2+2*n*n)){
        diagonal_counter++;
        theta[j]= arma::randn()*0.0001+theta_mu[j];
      }
      else{
        theta[j]=arma::randn()*0.00001+theta_mu[j];
      }
    }


    //arma::mat C = arma::zeros(n,n);
    int  index=0;
    for(int j=0; j <n; j++){
      for (int k = j;k < n; k++) {
        C(k,j)=theta[index];
        index++;
      }
    }

    A = arma::reshape(theta.rows((n * (n+1)/2), pow(n, 2) + (n * (n + 1)/2) - 1), n, n);
    B = arma::reshape(theta.rows(pow(n, 2) + (n * (n + 1)/2), 2*pow(n, 2) + (n * (n + 1)/2) - 1), n, n);
    G = arma::reshape(theta.rows(2*pow(n, 2) + (n * (n + 1)/2), 3*pow(n, 2) + (n * (n + 1)/2) - 1), n, n);

    if(valid_asymm_bekk(C,A,B,G,r,signs)){
      l++;
      double llv=loglike_asymm_bekk(theta,r,signs);

      if(llv>best_val){

        m++;
        best_val=llv;
        thetaOptim=theta;
        theta_mu=thetaOptim;

      }
      if(l>(2000/log(2+nc)) || m>=5){
        theta_mu=thetaOptim;
      }
    }

  }
  return Rcpp::List::create(Rcpp::Named("thetaOptim") = thetaOptim,
                            Rcpp::Named("best_val") = best_val);

}


// [[Rcpp::export]]
Rcpp::List sigma_bekk(arma::mat& r, arma::mat& C, arma::mat& A, arma::mat& G) {
// Computation of second order moment time paths and GARCH innovations
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat sigma = arma::zeros(r.n_rows, N2);
  arma::mat et = arma::zeros(r.n_rows, N);
  arma::mat ht = r.t() * r / r.n_rows;
  sigma.row(0) = arma::vectorise(ht).t();

  et.row(0) <- inv_gen(arma::sqrtmat_sympd(ht)) *  r.row(0).t();

  arma::mat CC  = C.t() * C;
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();

  for (int i = 1; i < r.n_rows; i++) {
    ht = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * ht * G;
    sigma.row(i) = arma::vectorise(ht).t();
    et.row(i) = (inv_gen(arma::chol(ht).t()) *  r.row(i).t()).t();
  }

  return Rcpp::List::create(Rcpp::Named("sigma_t")= sigma,
                            Rcpp::Named("e_t") = et);
}

// [[Rcpp::export]]
Rcpp::List sigma_bekk_asymm(arma::mat& r, arma::mat& C, arma::mat& A, arma::mat& B, arma::mat& G,arma::mat signs) {
  // Computation of second order moment time paths and GARCH innovations
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat sigma = arma::zeros(r.n_rows, N2);
  arma::mat et = arma::zeros(r.n_rows, N);
  arma::mat ht = r.t() * r / r.n_rows;
  sigma.row(0) = arma::vectorise(ht).t();

  et.row(0) <- inv_gen(arma::sqrtmat_sympd(ht)) * r.row(0).t();
//changed CC to C.t() * C instead of C * C.t() because C is upper triagular in the inputs
  arma::mat CC = C.t() * C;
  arma::mat At = A.t();
  arma::mat Bt = B.t();
  arma::mat Gt = G.t();

  for (int i = 1; i < r.n_rows; i++) {
    ht = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + indicatorFunction(r.row(i-1),signs)* Bt * r.row(i - 1).t() * r.row(i - 1) * B + Gt * ht * G;
    sigma.row(i) = arma::vectorise(ht).t();
    et.row(i) = (inv_gen(arma::chol(ht).t()) * r.row(i).t()).t();
  }

  return Rcpp::List::create(Rcpp::Named("sigma_t")= sigma,
                            Rcpp::Named("e_t") = et);
}

// [[Rcpp::export]]
arma::mat hesse_bekk(arma::mat theta, arma::mat r){
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
  arma::mat c0 = arma::zeros(N, N);
  int index =0;
  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      c0(j, i) = theta(index);
      index += 1;
    }
  }
  c0=c0.t();
  arma::mat a = arma::reshape(theta.rows((N * (N+1)/2), (pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat g = arma::reshape(theta.rows(((pow(N, 2) + (N * (N + 1)/2))), (2*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);

  arma::mat c_full = c0.t() * c0;



//often done calculations
      arma::mat t_kron_g=arma::kron(g,g).t();
      arma::mat gt=g.t();
      arma::mat at=a.t();
      arma::mat kron_comm_cgt=arma::kron(K_commutation,arma::reshape(gt,gt.size(),1));

      // Partial derivatives for initial period t = 1
      arma::mat ht = r.t() * r / r.n_rows;

      arma::mat dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), at * ht);
      arma::mat dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), gt * ht);
      arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N,N)) * L_elimination.t();

      arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

      arma::mat ht_inv = arma::inv(ht);
//Hessian
      arma::mat hessian = arma::zeros(theta.n_rows,theta.n_rows);

//Second derivatives for t=1
      arma::mat dHdada = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation));
      arma::mat dHdadc = arma::zeros(pow(N2,2),NoOfVars_C);
      arma::mat dHdadg = arma::zeros(pow(N2,2),N2);

      arma::mat dHdgda = arma::zeros(pow(N2,2),N2);
      arma::mat dHdgdg = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation));
      arma::mat dHdgdc = arma::zeros(pow(N2,2),NoOfVars_C);

      arma::mat dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t();
      arma::mat dHdcdg = arma::zeros(NoOfVars_C*N2,N2);
      arma::mat dHdcda = arma::zeros(NoOfVars_C*N2,N2);



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

              dHdada =C3* (arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(r.row(i-1).t()*r.row(i-1),arma::eye(N,N)))*K_commutation)+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdada;
              dHdadg = arma::kron(dHda.t(), arma::eye(N2,N2))*C1*(kron_comm_cgt + arma::kron(arma::reshape(gt,N2,1),K_commutation))+(arma::kron(arma::eye(N2,N2),arma::kron(g,gt).t()))*dHdadg;
              dHdadc = dHdadc; // Always zero (row may be deleted later on)

              dHdgda = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(arma::eye(N,N),gt)*dHda ))+(arma::kron(arma::eye(N2,N2),t_kron_g))*dHdgda;
              dHdgdg = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation+arma::kron(arma::eye(N,N),gt)*dHdg ))+arma::kron(dHdg.t(),arma::eye(N2,N2))*C1*(kron_comm_cgt+arma::kron(arma::reshape(gt,N2,1),K_commutation))+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdgdg;
              dHdgdc = C3*arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(arma::eye(N,N),gt)*dHdc)+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdgdc;

              dHdcda = dHdcda;// Always zero (row may be deleted later on)
              dHdcdg = arma::kron(dHdc.t(),arma::eye(N2,N2))*C1*(kron_comm_cgt+arma::kron(arma::reshape(gt,N2,1),K_commutation))+arma::kron(arma::eye(NoOfVars_C,NoOfVars_C),t_kron_g)*dHdcdg;
              dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t()+arma::kron(arma::eye(NoOfVars_C,NoOfVars_C),t_kron_g)*dHdcdc;

//updating first derivatives
                dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), a.t() * r.row(i-1).t()*r.row(i-1)) + arma::kron(g, g).t() * dHda;
                dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), gt * ht) + arma::kron(g, g).t() * dHdg;
                dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N,N)) * L_elimination.t() + arma::kron(g, g).t() * dHdc;

//update h
                ht = c_full + at * r.row(i-1).t()*r.row(i-1) * a + gt * ht *g;

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
arma::mat hesse_asymm_bekk(arma::mat theta, arma::mat r, arma::mat& signs){
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
  arma::mat c0 = arma::zeros(N, N);
  int index =0;
  for(int i = 0; i < N; i++){
    for (int j = i; j < N; j++) {
      c0(j, i) = theta(index);
      index += 1;
    }
  }
  c0=c0.t();
  arma::mat a = arma::reshape(theta.rows((N * (N+1)/2), (pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat b = arma::reshape(theta.rows(((pow(N, 2) + (N * (N + 1)/2))), (2*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat g = arma::reshape(theta.rows((( 2*pow(N, 2)  + (N * (N + 1)/2))), (3*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);

  arma::mat c_full = c0.t() * c0;



  //often done calculations
  arma::mat t_kron_g=arma::kron(g,g).t();
  arma::mat gt=g.t();
  arma::mat at=a.t();
  arma::mat bt=b.t();
  arma::mat kron_comm_cgt=arma::kron(K_commutation,arma::reshape(gt,gt.size(),1));

  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  arma::mat dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), at * ht);
  arma::mat dHdb = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), at * ht);
  arma::mat dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), gt * ht);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N,N)) * L_elimination.t();

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda,dHdb, dHdg).t();

  arma::mat ht_inv = arma::inv(ht);
  //Hessian
  arma::mat hessian = arma::zeros(theta.n_rows,theta.n_rows);

  //Second derivatives for t=1
  arma::mat dHdada = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation));
  arma::mat dHdadc = arma::zeros(pow(N2,2),NoOfVars_C);
  arma::mat dHdadg = arma::zeros(pow(N2,2),N2);
  arma::mat dHdadb = arma::zeros(pow(N2,2),N2);

  arma::mat dHdbdb = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation));
  arma::mat dHdbdc = arma::zeros(pow(N2,2),NoOfVars_C);
  arma::mat dHdbdg = arma::zeros(pow(N2,2),N2);
  arma::mat dHdbda = arma::zeros(pow(N2,2),N2);

  arma::mat dHdgda = arma::zeros(pow(N2,2),N2);
  arma::mat dHdgdb = arma::zeros(pow(N2,2),N2);
  arma::mat dHdgdg = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation));
  arma::mat dHdgdc = arma::zeros(pow(N2,2),NoOfVars_C);

  arma::mat dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t();
  arma::mat dHdcdg = arma::zeros(NoOfVars_C*N2,N2);
  arma::mat dHdcda = arma::zeros(NoOfVars_C*N2,N2);
  arma::mat dHdcdb = arma::zeros(NoOfVars_C*N2,N2);



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
  arma::mat GGt = arma::kron(g, g).t();
  for (int i=1; i<n;i++){

    dHdada =C3* (arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(r.row(i-1).t()*r.row(i-1),arma::eye(N,N)))*K_commutation)+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdada;
    dHdadg = arma::kron(dHda.t(), arma::eye(N2,N2))*C1*(kron_comm_cgt + arma::kron(arma::reshape(gt,N2,1),K_commutation))+(arma::kron(arma::eye(N2,N2),arma::kron(g,gt).t()))*dHdadg;
    dHdadc = dHdadc; // Always zero (row may be deleted later on)
    dHdadb = dHdadb; // Always zero (row may be deleted later on)

    dHdbdb = indicatorFunction(r.row(i-1),signs)*C3* (arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(r.row(i-1).t()*r.row(i-1),arma::eye(N,N)))*K_commutation)+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdbdb;
    dHdbdg = arma::kron(dHdb.t(), arma::eye(N2,N2))*C1*(kron_comm_cgt + arma::kron(arma::reshape(gt,N2,1),K_commutation))+(arma::kron(arma::eye(N2,N2),arma::kron(g,gt).t()))*dHdbdg;
    dHdbdc = dHdbdc; // Always zero (row may be deleted later on)
    dHdbda = dHdbda; // Always zero (row may be deleted later on)

    dHdgda = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(arma::eye(N,N),gt)*dHda ))+(arma::kron(arma::eye(N2,N2),t_kron_g))*dHdgda;
    dHdgdb = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(arma::eye(N,N),gt)*dHdb ))+(arma::kron(arma::eye(N2,N2),t_kron_g))*dHdgdb;
    dHdgdg = C3*(arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(ht,arma::eye(N,N))*K_commutation+arma::kron(arma::eye(N,N),gt)*dHdg ))+arma::kron(dHdg.t(),arma::eye(N2,N2))*C1*(kron_comm_cgt+arma::kron(arma::reshape(gt,N2,1),K_commutation))+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdgdg;
    dHdgdc = C3*arma::kron(arma::reshape(arma::eye(N,N),N2,1),arma::kron(arma::eye(N,N),gt)*dHdc)+arma::kron(arma::eye(N2,N2),t_kron_g)*dHdgdc;

    dHdcda = dHdcda;// Always zero (row may be deleted later on)
    dHdcdb = dHdcdb;// Always zero (row may be deleted later on)
    dHdcdg = arma::kron(dHdc.t(),arma::eye(N2,N2))*C1*(kron_comm_cgt+arma::kron(arma::reshape(gt,N2,1),K_commutation))+arma::kron(arma::eye(NoOfVars_C,NoOfVars_C),t_kron_g)*dHdcdg;
    dHdcdc = 2*arma::kron(L_elimination,D_duplication*D_gen_inv)*C1*arma::kron(arma::eye(N2,N2),arma::reshape(arma::eye(N,N),N2,1))*L_elimination.t()+arma::kron(arma::eye(NoOfVars_C,NoOfVars_C),t_kron_g)*dHdcdc;

    //updating first derivatives
    dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), a.t() * r.row(i-1).t()*r.row(i-1)) + arma::kron(g, g).t() * dHda;
    dHdb = indicatorFunction(r.row(i-1),signs) *  2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), bt * r.row(i-1).t() * r.row(i-1)) + GGt * dHdb;
    dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), gt * ht) + arma::kron(g, g).t() * dHdg;
    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N,N)) * L_elimination.t() + arma::kron(g, g).t() * dHdc;

    //update h
    ht = c_full + at * r.row(i-1).t()*r.row(i-1) * a + indicatorFunction(r.row(i-1),signs) *bt * r.row(i-1).t()*r.row(i-1) * b + gt * ht *g;

    //Computation of Hessian for each t>1


    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda,dHdb, dHdg).t();


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
arma::mat eigen_value_decomposition(arma::mat& A){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym( eigval, eigvec, A );


  arma::mat diag_mat_eigval = arma::diagmat(sqrt(eigval));
  return eigvec*diag_mat_eigval*eigvec.t();

}

/*
// [[Rcpp::export]]
Rcpp::List recursive_search_BEKK(arma::mat r, arma::vec c0, arma::vec avec,
                                arma::vec gvec, int index, arma::mat thetaopt, double likmax){

  int n  = r.n_cols;
  int start = -3;
  int endr = 3;
  int step = 6;
  int indextest = 0;

  if (index == pow(n, 2)){
    index += 2;
  } else if (index < pow(n, 2)){
    indextest = (index-1)/(n+1);
  } else{
    indextest = (index-pow(n, 2)-1)/(n+1);
  }
  // we have a diagonal element
  if (indextest - floor(indextest) == 0){
    index += 1;
  }

  arma::vec seq = arma::regspace(start, step, endr);

  for (int i = 0; i < seq.size(); i++){
    double inner = seq(i);
    double val = inner/100;

    // set a and g respectively according to index, exclude diagonal elements
    if (index <= pow(n, 2)) {
      avec(index-1) = val;
    } else {
      gvec(index-pow(n, 2)-1) = val;
    }
    // last element is excluded
    if (index < (2*pow(n, 2)-1)) {
      // recursive step
      Rcpp::List result = recursive_search_BEKK(r, c0, avec, gvec, index+1, thetaopt, likmax);
      arma::mat thetaopt = result(0);
      likmax = result(1);
    } else{
      // final step
      arma::mat theta = arma::join_cols(c0,avec,gvec);
      //likelihood
      double lik = loglike_bekk(theta, r);
      if (lik > likmax){
        thetaopt = theta;
        likmax = lik;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("theta")= thetaopt,
                              Rcpp::Named("lik") = likmax);
}

// [[Rcpp::export]]
Rcpp::List grid_search_BEKK(arma::mat r) {
  int N  = r.n_cols;

  arma::mat uncond_var = r.t() * r / r.n_rows;
  arma::mat A = arma::zeros(N, N);
  arma::mat G = arma::zeros(N, N);
  arma::mat C = arma::zeros(N, N);

  arma::mat hel = 0.05*uncond_var.diag();


  for (int i = 0; i < N; i++) {
    A(i,i) = 0.3;
    G(i,i) = 0.92;
    C(i,i) = hel(i,0);
  }

  double cij;

  for (int i = 0; i < N; i++){
    for (int j = i; j < N; j++){
      cij = uncond_var(i, j)/sqrt(uncond_var(i,i)*uncond_var(j,j));
      C(i,j) = cij*sqrt(C(i, i)*C(j, j));
      C(j, i) = C(i,j);
    }
  }

  C = chol(C, "lower");

  arma::mat C0;
  C0 = C.col(0);

  if (N > 2) {
    for (int i = 1; 1 < (N-1); i++){
      C0 = arma::join_cols(C0, C(arma::span(i,N-1), i));
    }
  }

  int csize = C0.size();
  C0.resize(csize+1);
  C0(csize) = C(N-1,N-1);


  // deterministic variance components

  arma::mat th0 = arma::join_cols(C0, arma::vectorise(A), arma::vectorise(G));

  // change elements of A and G and compute likelihood in each step

  double likmax = -1e25;

  Rcpp::List result = recursive_search_BEKK(r, C0, arma::vectorise(A), arma::vectorise(G), 1, th0, likmax);
  //th0 = result[0];
  //likmax=result[[2]]
  //return(list(th0,likmax))

  return result;
}
 */
