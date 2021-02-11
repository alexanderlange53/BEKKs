#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::mat elimination_mat(int n) {
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
arma::mat commutation_mat(int n) {
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
arma::mat duplication_mat(int n) {
  // Generates a duplication matrix for size 'n'
  int n2 = pow(n, 2);

  arma::mat el = elimination_mat(n);
  arma::mat co = commutation_mat(n);
  arma::mat m = arma::eye(n2, n2) + co;

  arma::mat dup = m*el.t()*arma::inv(el*m*el.t());

  return dup;
}

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

// [[Rcpp::export]]
arma::mat score_bekk(arma::mat theta, arma::mat r) {
  int N = r.n_cols;
  int N2 = pow(N, 2);

  arma::mat L_elimination = elimination_mat(N);
  arma::mat D_duplication = duplication_mat(N);


  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();

  arma::mat gradients = arma::zeros(r.n_rows, theta.n_rows);


  // transform vector of parameter 'theta' to actual parameter matrices of BEKK model
  arma::mat c0 = arma::zeros(N, N);
  arma::uvec lw_idx = arma::trimatu_ind(arma::size(c0));
  c0.elem(lw_idx) = theta.rows(0, (N * (N + 1)/2) - 1);
  arma::mat a = arma::reshape(theta.rows((N * (N+1)/2), (pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);
  arma::mat g = arma::reshape(theta.rows(((pow(N, 2) + (N * (N + 1)/2))), (2*pow(N, 2) + (N * (N + 1)/2) - 1)), N, N);

  arma::mat c_full = c0.t() * c0;


  // Partial derivatives for initial period t = 1
  arma::mat ht = r.t() * r / r.n_rows;

  arma::mat dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), a.t() * ht);
  arma::mat dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N,N), g.t() * ht);
  arma::mat dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N,N)) * L_elimination.t();

  arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

  arma::mat ht_sqrt_inv = arma::inv(arma::sqrtmat_sympd(ht));

  arma::vec et = ht_sqrt_inv * r.row(0).t();
  // Score function
  for (int k = 0; k < theta.n_rows; k++) {
    arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);

    arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());

    gradients(0, k) = -(0.5) * arma::sum(mat_temp.diag());
  }

  // Partial derivatives for period t >= 2
  arma::mat GGt = arma::kron(g, g).t();
  for (int i = 1; i < r.n_rows; i++) {
    dHda = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), a.t() * r.row(i-1).t() * r.row(i-1)) + GGt * dHda;
    dHdg = 2 * D_duplication * D_gen_inv * arma::kron(arma::eye(N, N), g.t() * ht) + GGt * dHdg;
    dHdc = 2 * D_duplication * D_gen_inv * arma::kron(c0.t(), arma::eye(N, N)) * L_elimination.t() + GGt * dHdc;

    arma::mat dHdtheta = arma::join_horiz(dHdc, dHda, dHdg).t();

    ht = c_full + a.t() * r.row(i-1).t() * r.row(i-1) * a + g.t() * ht * g;

    ht_sqrt_inv = arma::inv(arma::sqrtmat_sympd(ht));
    et = ht_sqrt_inv * r.row(i).t();

    for (int k = 0; k < theta.n_rows; k++) {
      arma::mat dh = arma::reshape(dHdtheta.row(k), N, N);

      arma::mat mat_temp = ht_sqrt_inv * dh * ht_sqrt_inv * (arma::eye(N, N) - et * et.t());

      gradients(i, k) = -(0.5) * arma::sum(mat_temp.diag());
    }
  }

  return gradients;
}
