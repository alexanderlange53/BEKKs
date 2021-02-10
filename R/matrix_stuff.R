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

  return(dup)
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
