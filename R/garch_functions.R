# Generates vector of GRACH variances (Sigma_e)
greatZ0_garch <- function(Tob, Z, q, p, ucvar, theta) {
  theta <- matrix(theta, nrow = p+q+1, ncol = 1)
  vvec <- rep(0, Tob-q)

  e <- ucvar * rep(1, p)
  vvec <- c(e, vvec)

  g <- nrow(theta)

  for (i in (p+1):(Tob+p)) {
    b <- vvec[(i-p):(i-1)]
    b <- rev(b)
    b <- theta[(q+2):g]*b
    b <- sum(b)

    vvec[i] <- Z[i-p, ] %*% theta[1:(q+1)] + b
  }

  return(vvec[(p+1):(Tob+p)])
}
