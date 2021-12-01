process_object <- function(x) {
  UseMethod('process_object')
}

process_object.bekkFit <- function(x) {
  theta <- x$theta
  N <- ncol(x$C0)

  return(list(theta = theta,
              N = N))
}

process_object.bekkSpec <- function(x) {
  if (is.null(x$init_values)) {
    stop('Please provide "initial_values" in "bekk_spec" as paramater for simulation.')
  }
  if (is.null(x$N)) {
    stop('Please provide "N" in "bekk_spec" as dimension for simulation.')
  }

  theta <- x$init_values
  N <- x$N

  return(list(theta = theta,
              N = N))
}

# Obtaining QML t-ratios
QML_t_ratios <- function(theta, r) {
  s1 <- score_bekk(theta, r)
  s1 <- crossprod(s1)

  s2 <- hesse_bekk(theta, r)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(theta/s2)
}

QML_t_ratios_asymm <- function(theta, r, signs) {
  s1 <- score_asymm_bekk(theta, r, signs)
  s1 <- crossprod(s1)

  s2 <- hesse_asymm_bekk(theta, r, signs)
  s2 <- solve(s2) %% s1 %% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(theta/s2)
}


simulate_bekk_alternative=function(theta,NoOBs,n){

  #   #Length of each series
  series=matrix(0,ncol = n,nrow=NoOBs)
  numb_of_vars <- 2*n^2+n*(n+1)/2
  C <- matrix(0,ncol = n,nrow = n)
  index <- 1
  for(i in 1:n){
    for (j in i:n) {
      C[j,i] <- theta[index]
      index <- index+1
    }
  }
  C <- C
  C_full=crossprod(C)
  A = matrix(theta[index:(index+n^2-1)], n)
  At=t(A)
  G = matrix(theta[(index+n^2):numb_of_vars], n)
  Gt=t(G)
  cov_mat=cor(BiCopSim(100000,4,1.3))
  cov_dec=solve(t(chol(cov_mat)))

  #unconditional variance
  Uncond_var=matrix(solve(diag(n^2) - t(kronecker(A, A)) - t(kronecker(G, G))) %*% c(C_full),n)

  H=Uncond_var
  h_dec=t(chol(H))
  series[1,]=h_dec%*%cov_dec%*%t(qnorm(BiCopSim(1,4,1.3)))
  for(i in 2:NoOBs){
    H=C_full+At%*%series[i-1,]%*%t(series[i-1,])%*%A+Gt%*%H%*%G
    h_dec=t(chol(H))
    series[i,]=h_dec%*%cov_dec%*%t(qnorm(BiCopSim(1,4,1.3)))
  }

  return(series)
}



