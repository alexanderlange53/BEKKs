bhh_bekk <- function(r, theta, max_iter) {

  steps <- seq(5,0, by = -0.5)
  step <- 0.01
  count_loop <- 1
  theta_candidate <- theta


  while (count_loop < max_iter) {
    theta <- theta_candidate
    theta_temp <- matrix(0, nrow = nrow(theta), length(steps))

    score_functio <- score_bekk(r)
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
  g <- c_mat$b

  c_full <- crossprod(c0)

  # Partial derivatives for initial perdio t = 1
  ht <- crossprod(r)/nrow(r)

  dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% h)
  dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% h)
  dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(c0, diag(N)) %*% t(L_elimination)

  dHdtheta <- cbind(dHda, dHdg, dHdc)

  ht_sqrt_inv <- solve(sqrtm(ht))

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
  a_0 <- theta[(N * (N+1)/2 + 1):(3 * N^2 + N)/2, ]
  g_0 <- theta[((3 * N^2 + N)/2 + 1):(2*N^2 + (N * (N + 1)/2)), ]

  a <- t(matrix(c_0, n, n))
  g <- t(matrix(g_0, n, n))

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              g = g))
}
