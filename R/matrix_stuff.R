# Converts parameter vector into coefficients matrices of BEKK model
coef_mat <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a_0 <- theta[(N * (N+1)/2 + 1):(N^2 + (N * (N + 1)/2)), ]
  g_0 <- theta[((N^2 + (N * (N + 1)/2)) + 1):(2*N^2 + (N * (N + 1)/2)), ]

  a <- (matrix(a_0, N, N))
  g <- (matrix(g_0, N, N))

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              g = g))
}
coef_mat_asymm <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a_0 <- theta[(N * (N+1)/2 + 1):(N^2 + (N * (N + 1)/2)), ]
  b_0 <- theta[(N^2 + (N * (N + 1)/2) + 1):(2*N^2 + (N * (N + 1)/2)), ]
  g_0 <- theta[((2*N^2 + (N * (N + 1)/2)) + 1):(3*N^2 + (N * (N + 1)/2)), ]

  a <- (matrix(a_0, N, N))
  b <- (matrix(b_0, N, N))
  g <- (matrix(g_0, N, N))

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              b = b,
              g = g))
}

coef_mat_diagonal <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a_0 <- theta[(N * (N+1)/2 + 1):(N + (N * (N + 1)/2)), ]
  g_0 <- theta[((N + (N * (N + 1)/2)) + 1):(2*N + (N * (N + 1)/2)), ]

  a <- diag(a_0)
  g <- diag(g_0)

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              g = g))
}
coef_mat_asymm_diagonal <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a_0 <- theta[(N * (N+1)/2 + 1):(N + (N * (N + 1)/2)), ]
  b_0 <- theta[(N + (N * (N + 1)/2) + 1):(2*N + (N * (N + 1)/2)), ]
  g_0 <- theta[((2*N + (N * (N + 1)/2)) + 1):(3*N + (N * (N + 1)/2)), ]

  a <- diag(a_0)
  b <- diag(b_0)
  g <- diag(g_0)

  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              b = b,
              g = g))
}

coef_mat_scalar <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a <- theta[(N * (N+1)/2 + 1), ]
  g <- theta[(( (N * (N + 1)/2)) + 2), ]



  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              g = g))
}
coef_mat_asymm_scalar <- function(theta, N) {
  c_0 <- theta[1:(N * (N + 1)/2), ]
  a <- theta[(N * (N+1)/2 + 1), ]
  b <- theta[(N * (N+1)/2 + 2), ]
  g <- theta[(( (N * (N + 1)/2)) + 3), ]



  c0 <- matrix(0, N, N)
  c0[lower.tri(c0, diag = TRUE)] <- c_0

  return(list(c0 = c0,
              a = a,
              b = b,
              g = g))
}
