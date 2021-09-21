gridSearch_BEKK <- function(r) {
  N  <- ncol(r)

  uncond_var <- crossprod(r)/nrow(r)
  A <- matrix(0, ncol = N, nrow = N)
  G <- matrix(0, ncol = N, nrow = N)
  C <- matrix(0, ncol = N, nrow = N)

  diag(A) <- 0.3
  diag(G) <- 0.92
  diag(C) <- 0.05*diag(uncond_var)


  for (i in 1:N){
    for (j in seq(i,N)){

      cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
      C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
      C[j,i] <- C[i, j]

    }
  }

  C = t(chol(C))
  C0 = C[,1]

  if (N > 2) {
    for (i in 2:(N-1)){
      C0 = c(C0, C[i:N, i])
    }
  }

  C0 = c(C0, C[N, N])

  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j) {
        A[i, j] <-  0.03
        if (j == N & i == 1) {
          A[i, j] <-  -0.03
        }
        G[i, j] <- -0.03
      } else if (i > j) {
        A[i, j] <- 0.03
        G[i, j] <- 0.03
      }
    }
  }

  th0 = c(C0, c(A), c(G))
  lik = loglike_bekk(th0, r)

  return(list(th0, lik))

}

gridSearch_asymmetricBEKK <- function(r) {
  N  <- ncol(r)

  uncond_var <- crossprod(r)/nrow(r)
  A <- matrix(0, ncol = N, nrow = N)
  B <- matrix(0, ncol = N, nrow = N)
  G <- matrix(0, ncol = N, nrow = N)
  C <- matrix(0, ncol = N, nrow = N)
  #th0=numeric(2*n^2+n*(n+1)/2)

  diag(A) <- 0.25
  diag(B) <- 0.05
  diag(G) <- 0.92
  diag(C) <- 0.05*diag(uncond_var)


  for (i in 1:N){
    for (j in seq(i,N)){

      cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
      C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
      C[j,i] <- C[i, j]

    }
  }

  C = t(chol(C))
  C0 = C[,1]

  if (N > 2) {
    for (i in 2:(N-1)){
      C0 = c(C0, C[i:N, i])
    }
  }

  C0 = c(C0, C[N, N])

  for (i in 1:N) {
    for (j in 1:N) {
      if (i < j) {
        A[i, j] <-  0.02
        B[i, j] <-  0.01
        if (j == N & i == 1) {
          A[i, j] <-  -0.02
          B[i, j] <-  -0.01
        }
        G[i, j] <- -0.03
      } else if (i > j) {
        A[i, j] <- 0.02
        B[i, j] <- 0.01
        G[i, j] <- 0.03
      }
    }
  }

  th0 = c(C0, c(A), c(B), c(G))
  lik = loglike_asymm_bekk(th0, r)

  return(list(th0, lik))

}
