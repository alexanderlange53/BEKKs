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

gridSearch_asymmetricBEKK <- function(r, signs) {
  N  <- ncol(r)

  uncond_var <- crossprod(r)/nrow(r)
  A <- matrix(0, ncol = N, nrow = N)
  B <- matrix(0, ncol = N, nrow = N)
  G <- matrix(0, ncol = N, nrow = N)
  C <- matrix(0, ncol = N, nrow = N)
  #th0=numeric(2*n^2+n*(n+1)/2)

  diag(A) <- 0.3
  diag(B) <- 0.3
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
        B[i, j] <-  0.03
        if (j == N & i == 1) {
          A[i, j] <-  -0.03
          B[i, j] <-  -0.03
        }
        G[i, j] <- -0.03
      } else if (i > j) {
        A[i, j] <- 0.03
        B[i, j] <- 0.03
        G[i, j] <- 0.03
      }
    }
  }

  th0 = c(C0, c(A), c(B), c(G))
  lik = loglike_asymm_bekk(th0, r, signs)

  return(list(th0, lik))

}

gridSearch_dBEKK <- function(r) {
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

  th0 = c(C0, diag(A), diag(G))
  lik = loglike_dbekk(th0, r)

  return(list(th0, lik))

}

gridSearch_asymmetricdBEKK <- function(r, signs) {
  N  <- ncol(r)

  uncond_var <- crossprod(r)/nrow(r)
  A <- matrix(0, ncol = N, nrow = N)
  B <- matrix(0, ncol = N, nrow = N)
  G <- matrix(0, ncol = N, nrow = N)
  C <- matrix(0, ncol = N, nrow = N)
  #th0=numeric(2*n^2+n*(n+1)/2)

  diag(A) <- 0.3
  diag(B) <- 0.3
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


  th0 = c(C0, diag(A), diag(B), diag(G))
  lik = loglike_asymm_dbekk(th0, r, signs)

  return(list(th0, lik))

}

gridSearch_sBEKK <- function(r) {
  N  <- ncol(r)
  C <- matrix(0, ncol = N, nrow = N)
  uncond_var <- crossprod(r)/nrow(r)


  a <- 0.2
  g <- 0.7
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



  th0 = c(C0, a, g)
  lik = loglike_sbekk(th0, r)

  return(list(th0, lik))

}

gridSearch_asymmetricsBEKK <- function(r, signs) {
  N  <- ncol(r)
  C <- matrix(0, ncol = N, nrow = N)
  uncond_var <- crossprod(r)/nrow(r)

  #th0=numeric(2*n^2+n*(n+1)/2)

  a <- 0.2
  b <- 0.1
  g <- 0.6
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



  th0 = c(C0, a, b, g)
  lik = loglike_asymm_sbekk(th0, r, signs)

  return(list(th0, lik))

}
#H2 Code
# gridSearch_BEKK <- function(r){
#   N  <- ncol(r)
#
#   uncond_var <- crossprod(r)/nrow(r)
#   A <- matrix(0, ncol = N, nrow = N)
#   G <- matrix(0, ncol = N, nrow = N)
#   C <- matrix(0, ncol = N, nrow = N)
#   #th0=numeric(2*n^2+n*(n+1)/2)
#
#   diag(A) <- 0.3
#   diag(G) <- 0.92
#   diag(C) <- 0.05*diag(uncond_var)
#
#
#   for (i in 1:N){
#     for (j in seq(i,N)){
#
#       cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
#       C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
#       C[j,i] <- C[i, j]
#
#     }
#   }
#
#   C = t(chol(C))
#   C0 = C[,1]
#
#   if (N > 2) {
#     for (i in 2:(N-1)){
#       C0 = c(C0, C[i:N, i])
#     }
#   }
#
#   C0 = c(C0, C[N, N])
#
#   # deterministic variance components
#
#   th0 = c(C0, c(A), c(G))
#
#   # change elements of A and G and compute likelihood in each step
#
#   likmax = -1e25
#
#   #print th0
#   result= recursiveSearch_BEKK(r, C0, c(A), c(G), 1, th0, likmax)
#   th0=result[[1]]
#   likmax=result[[2]]
#   return(list(th0,likmax))
#
# }
#
#
# recursiveSearch_BEKK=function(r, c0, avec, gvec, index, thetaopt, likmax){
#
#   n=ncol(r)
#   start= -3
#   endr = 3
#   step = 6
#   indextest=0
#
#   if (index == n^2){
#     index = index + 2;
#   } else if (index < n^2){
#     indextest = (index-1)/(n+1)
#   } else{
#     indextest = (index-n^2-1)/(n+1)
#   }
#   # we have a diagonal element
#   if (indextest - floor(indextest) == 0){
#     index = index + 1
#   }
#
#   for (i in seq(start, endr, step)){
#     val = i/100
#
#     #set a and g respectively according to index, exclude diagonal elements
#     if (index <= n^2){
#       avec[index] = val
#     } else{
#       gvec[index-n^2] = val
#     }
#     #last element is excluded
#     if (index < (2*n^2-1)){
#       # recursive step
#       result= recursiveSearch_BEKK(r, c0, avec, gvec, index+1, thetaopt, likmax)
#       thetaopt=result[[1]]
#       likmax=result[[2]]
#     } else{
#       #final step
#       theta = c(c0,avec,gvec)
#       #likelihood
#       lik = loglike_bekk(theta,r)
#       if (lik > likmax){
#         thetaopt = theta
#         likmax = lik
#       }
#
#     }
#   }
#
#   return(list(thetaopt,likmax))
# }

