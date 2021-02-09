bhh_bekk <- function(r, theta, max_iter) {

  steps <- seq(5,0, by = -0.5)
  step <- 0.01
  count_loop <- 1
  theta_candidate <- theta


  while (count_loop < max_iter) {
    theta <- theta_candidate
    theta_temp <- matrix(0, nrow = nrow(theta), length(steps))

    score_function <- score_bekk(theta, r)
    outer_score <- tcrossprod(score_function)
    score_function <- colSums(score_function)

    # Hier likelihood
  }


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

# Generates an elimination matrix for size 'n'
elimination_mat <- function(n) {
  n1 <- n * (n + 1) / 2
  n2 <- n^2

  init <- diag(n1)
  oes <- 1

  eli <- matrix(init[, 1], ncol = 1)

  block <- n

  while (ncol(eli) < n2) {
    if (ncol(eli) == 1) {
      eli <- init[, 1:block]
    } else {
      eli <- cbind(eli, init[, 1:block])
    }

    if (ncol(init) > 1) {
      init <- matrix(init[, (block+1):ncol(init)], nrow = n1)
    }

    eli <- cbind(eli, (matrix(0, nrow = nrow(eli), ncol = oes)))

    oes <- oes + 1

    block <- block - 1
  }

  eli <- eli[, 1:(ncol(eli)-n)]
  return(eli)
}

# Generates a duplication matrix for size 'n'
duplication_mat <- function(n) {
  n2 <- n^2

  el <- elimination_mat(n)
  co <- commutation_mat(n)
  m <- diag(n2) + co

  dup <- m%*%t(el)%*%solve(el%*%m%*%t(el))

}

# generates a (square) commutation matrix for 'n'
commutation_mat <- function(n) {
  K <- matrix(0, n^2, n^2)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i + n*(j - 1), j + n*(i - 1)] <- 1
    }
  }
  return(K)
}

# Converts parameter vefctor into coefficients matrices of BEKK model
coef_mat <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a_0 <- theta[(N * (N+1)/2 + 1):(N^2 + (N * (N + 1)/2)), ]
  g_0 <- theta[((N^2 + (N * (N + 1)/2)) + 1):(2*N^2 + (N * (N + 1)/2)), ]

  a <- (matrix(a_0, N, N))
  g <- (matrix(g_0, N, N))

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = t(c0),
              a = a,
              g = g))
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

#Computation of BEKK
comph_bekk <- function(C, A, G, r) {
  # dimensions
  n      = nrow(r)

  # redefine for convenience
  CC     = t(C) %*% C
  At     = t(A)
  Gt     = t(G)

  # compute uncond. covariance matrix
  Hu     = (t(r) %*% r) / nrow(r)

  # compute H trajectory
  H      = vector(mode = "list", n)
  H[[1]] = Hu
  for (i in 2:n) {
    H[[i]] <- (CC + At %*% (r[i - 1,] %*% t(r[i - 1,])) %*% A
               + Gt %*% H[[i - 1]] %*% G)
  }
  return(H)
}
#Log-Likelihood Value
loglike_bekk <- function(theta, r) {
  # convert to matrices
  n=ncol(r)
  #Length of each series
  NoOBs=nrow(r)
  numb_of_vars=2*n^2+n*(n+1)/2
  C=matrix(0,ncol = n,nrow = n)
  index=1
  for(i in 1 : n){
    for (j in i:n) {
      C[j,i]=theta[index]
      index=index+1
    }
  }
  A = matrix(theta[index:(index+n^2-1)], n)
  G = matrix(theta[(index+n^2):numb_of_vars], n)

  # check constraints
  if (!valid_bekk(C, A, G)) {
    return(-10*n)
  }

  # compute H
  H   = comph_bekk(C, A, G, r)
  # compute llv

  llv = numeric(NoOBs)
  for (i in 1:NoOBs) {
    llv[i] = c(t(r[i,]) %*% solve(H[[i]]) %*% r[i,])
  }
  llv = -0.5 * (2 * log(2 * pi) + log(sapply(H, det)) + llv)
  return(sum(llv))
}
