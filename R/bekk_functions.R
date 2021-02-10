bhh_bekk <- function(r, theta, max_iter, crit) {

  steps <- c(5,2,1,0.5,0.25,0.1,0.01,0.005,0)
  step <- 0.01
  count_loop <- 1
  theta_candidate <- theta
  exit_loop <- 0


  while (count_loop < max_iter & exit_loop == 0) {
    theta <- theta_candidate
    theta_temp <- matrix(0, nrow = nrow(theta), length(steps))

    score_function <- score_bekk(theta, r)
    outer_score <- crossprod(score_function)
    score_function <- colSums(score_function)

    # Hier likelihood
    lik <- loglike_bekk(theta, r)

    for (i in 1:length(steps)) {
      temp <- theta_candidate + step * steps[i]*solve(outer_score) %*% score_function # hier nochmal gucken H2 hat da kontrolliert on inverser exstiert
      theta_temp[, i] <- temp
    }

    likelihood_candidates <- rep(0, length(steps))
    likelihood_candidates[length(steps)] <- lik

    j <- length(steps) - 1
    exit_inner <-0
    while (j >= 1 & exit_inner == 0) {
      likelihood_candidates[j] <- loglike_bekk(theta_temp[, j], r)
      if (likelihood_candidates[j+1] > likelihood_candidates[j]) {
        exit_inner <-1
      }
      j <- j-1
    }

    max_index <- which.max(likelihood_candidates[(j+1):length(likelihood_candidates)]) + j
    likelihood_best <- likelihood_candidates[max_index]

    # exit criterion strange
    if ((likelihood_best - likelihood_candidates[length(likelihood_candidates)])^2/abs(likelihood_candidates[length(likelihood_candidates)]) < crit) {
      exit_loop <- 1
    }
    theta_candidate <- matrix(theta_temp[, max_index], ncol= 1)
    count_loop <- count_loop + 1
  }

  likelihood_final <- loglike_bekk(theta_candidate, r)
  score_final <- score_bekk(theta_candidate, r)
  s1 <- sqrt(diag(solve(crossprod(score_final))))

  t_val <- theta_candidate/s1

  return(list(theta = theta_candidate,
              t_val = t_val,
              likelihood = likelihood_final,
              iter = count_loop))
}

score_bekk <- function(theta, r) {
  N <- ncol(r)
  N2 <- N^2

  L_elimination <- elimination_mat(N)
  D_duplication <- duplication_mat(N)

  D_gen_inv <- solve(crossprod(D_duplication)) %*% t(D_duplication)

  gradients <- matrix(0, nrow = nrow(r), ncol = nrow(theta))

  c_mat <-  coef_mat(theta, N)

  c0 <- c_mat$c0
  a <- c_mat$a
  g <- c_mat$g

  c_full <- crossprod(c0)

  # Partial derivatives for initial period t = 1
  ht <- crossprod(r)/nrow(r)

  dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% ht)
  dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% ht)
  dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination)

  dHdtheta <- t(cbind(dHdc, dHda, dHdg))

  ht_sqrt_inv <- solve(sqrtm(ht))

  et <- ht_sqrt_inv %*% r[1,]
  # Score function
  for (k in 1:nrow(theta)) {
    dh <- matrix(dHdtheta[k, ], N, N)

    mat_temp <- ht_sqrt_inv %*% dh %*% ht_sqrt_inv %*% (diag(N) - tcrossprod(et))

    gradients[1, k] <- -(0.5) * sum(diag(mat_temp))
  }

  # Partial derivatives for period t >= 2
  for (i in 2:nrow(r)) {
    dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% r[(i-1), ] %*% t(r[(i-1), ])) + t(kronecker(g, g)) %*% dHda
    dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% ht) + t(kronecker(g, g)) %*% dHdg
    dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination) + t(kronecker(g, g)) %*% dHdc

    dHdtheta <- t(cbind(dHdc, dHda, dHdg))

    ht <- c_full + t(a) %*% r[(i-1), ] %*% t(r[(i-1), ]) %*% a + t(g) %*% ht %*%g

    ht_sqrt_inv <- solve(sqrtm(ht))
    et <- ht_sqrt_inv %*% r[i,]

    for (k in 1:nrow(theta)) {
      dh <- matrix(dHdtheta[k, ], N, N)

      mat_temp <- ht_sqrt_inv %*% dh %*% ht_sqrt_inv %*% (diag(N) - tcrossprod(et))

      gradients[i, k] <- -(0.5) * sum(diag(mat_temp))
    }
  }

  return(gradients)
}


#Testing existance, uniqueness and stationarity
valid_bekk <- function(C, A, G) {
  # condition for positive-definit covariance
  if (any(diag(C) < 0)) {
    return(FALSE)
  }

  # condition for uniqueness
  if (A[1, 1] < 0 || G[1, 1] < 0) {
    return(FALSE)
  }

  # check stationarity for BEKK(1,1): Engle & Kroner (1995), Prop. 2.7
  if (!all(abs(eigen(kronecker(A, A)
                     + kronecker(G, G))$values) < 1)) {
    return(FALSE)
  }

  return(TRUE)
}


# Log-Likelihood function
loglike_bekk <- function(theta, r) {
  # convert to matrices
  n <- ncol(r)
  #Length of each series
  NoOBs <- nrow(r)
  numb_of_vars <- 2*n^2+n*(n+1)/2
  C <- matrix(0,ncol = n,nrow = n)
  index <- 1
  for(i in 1:n){
    for (j in i:n) {
      C[j,i] <- theta[index]
      index <- index+1
    }
  }
  C <- t(C)
  A = matrix(theta[index:(index+n^2-1)], n)
  G = matrix(theta[(index+n^2):numb_of_vars], n)

  # check constraints
  if (!valid_bekk(C, A, G)) {
    return(-1e25)
  }

  # compute inital H
  H <- (t(r) %*% r) / nrow(r)

  CC  <- t(C) %*% C
  At  <- t(A)
  Gt  <- t(G)


  llv = numeric(NoOBs)
  llv[1] <- -0.5 * (n * log(2 * pi) + log(det(H)) + c(t(r[1,]) %*% solve(H) %*% r[1,]))
  for (i in 2:NoOBs) {
    H <- (CC + At %*% (r[i - 1,] %*% t(r[i - 1,])) %*% A + Gt %*% H %*% G)
    llv[i] = -0.5 * (n * log(2 * pi) + log(det(H)) + c(t(r[i,]) %*% solve(H) %*% r[i,]))
  }
  return(sum(llv))
}
