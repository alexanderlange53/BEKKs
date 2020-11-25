# Script containing the functions to estimate the GARCH(p, q) parameter

BHHH_garch <- function(r2, q, p, theta, epsilon2, Z, Tob, max_iter, crit, ucvar){
  # BHHH algorithm based on first order derivatives and outer product for ARCH(q) model

  step_width <- seq(1, 0, length.out = 100) # find the optimal step width

  quit <- 0
  iter <- 1

  while (quit == 0 & iter < max_iter) {

    score_function <- score_garch(epsilon2, Z, Tob, q, p, theta, ucvar)

    out_prod_grad <- score_function %*% t(score_function)

    sum_grad <- rowSums(score_function)

    theta_hat <- matrix(0, nrow = nrow(theta), ncol = 100)
    for (i in 1:100) {
      theta_hat[,i] <- theta + step_width[i] * solve(out_prod_grad) %*% sum_grad
    }

    theta_hat <- abs(theta_hat)

    candidates <- matrix(0, 100, 1)

    candidates[1, ] <- likelihood_garch(Z, Tob, q, p, theta, epsilon2, ucvar)

    quit_inner <- 0
    iter_inner <- 2

    while (iter_inner < 100 & quit_inner == 0) {
      candidates[iter_inner, ] <- likelihood_garch(Z, Tob, q, p, theta_hat[, iter_inner], epsilon2, ucvar)

      # if (candidates[iter_inner] < candidates[iter_inner-1]) {
      #   quit_inner <- 1
      #   cat('stop')
      # }
      iter_inner <- iter_inner + 1
    }

    min_lik <- which.min(candidates)

    best_candidate <- candidates[min_lik,]

    if ((best_candidate - candidates[1,])^2 / abs(candidates[1,]) <  crit) {
      quit <- 1
    } else {
      theta <- as.matrix(theta_hat[, min_lik])
      iter <- iter + 1
    }
  }

  lik_optim <- likelihood_garch(Z, Tob, q, p, theta, epsilon2, ucvar)

  return(c(theta, lik_optim))
}

score_garch <- function(epsilon2, Z, Tob, q, p, theta, ucvar) {
  theta <- matrix(theta, nrow = p+q+1, ncol = 1)
  vvec <- greatZ_garch(Tob, Z, q, p, ucvar, theta) # GARCH variance sigma^2
  vvec2 <- c(ucvar, ucvar, vvec) # combining with unconditonal variance
  vvec2 <- vvec2[1:(Tob-q)] # discarding last q observations
  vvec1 <- c(ucvar, vvec)
  vvec1 <- vvec1[1:(Tob-q)]

  aa <- ucvar/(1 - theta[(q + 2):nrow(theta)] %*% rep(1, nrow(theta) - q -q))
  aa2 <- 1/(1 - theta[(q+2):nrow(theta)] %*% rep(1, nrow(theta) - q -q))

  gz <- cbind(Z, vvec1, vvec2)
  nparam <- nrow(theta)

  gz <- gz[,1:nparam]
  vabl <- matrix(0, nrow=(Tob-q), ncol = nparam)

  e <- c(aa) * matrix(1, p, nparam)
  e[,1] <- c(aa2) * matrix(1, p, 1)
  vabl <- rbind(e, vabl)

  for (i in (p+1):(Tob-q+p)) {
    b <- matrix(vabl[(i - p):(i-1),], nrow = p)
    b <- b[rev(1:nrow(b)),]
    b <- theta[(q+2):nparam,]*b
    if(p > 1) {
      b <- t(colSums(b))
    }else{
      b <- t(b)
    }

    vabl[i,] <- gz[(i-p),] + b
  }

  vabl <- vabl[(p+1):(Tob-q+p),]

  siggi1 <- vabl/vvec
  siggi2 <- vvec^2
  siggi2 <- epsilon2/siggi2
  siggi2 <- matrix(rep(siggi2, ncol(vabl)), nrow= length(siggi2), ncol=ncol(vabl), byrow = FALSE)*vabl

  ltv <- 0.5 * (siggi2 - siggi1)

  return(t(ltv))
}

likelihood_garch <- function(Z, Tob, q, p, theta, epsilon2, ucvar) {
  # likelihood function for specific parameter 'theta'

  theta <- matrix(theta, nrow = p+q+1, ncol = 1)

  #sigma2 <- Z %*% theta
  sigma2 <- greatZ_garch(Tob, Z, q, p, ucvar, theta)

  term2 <- epsilon2 / sigma2

  loglik <- 0.5 * (Tob - q) * log(2 * pi) + 0.5 * sum(log(sigma2)) + 0.5 * sum(term2)
  return(loglik)
}


# Generates vector of GRACH variances (Sigma_e)
greatZ_garch <- function(Tob, Z, q, p, ucvar, theta) {
  vvec <- rep(0, Tob-q)

  e <- ucvar * rep(1, p)
  vvec <- c(e, vvec)

  g <- nrow(theta)

  for (i in (p+1):(Tob-q+p)) {
    b <- vvec[(i-p):(i-1)]
    b <- rev(b)
    b <- theta[(q+2):g]*b
    b <- sum(b)

    vvec[i] <- Z[i-p, ] %*% theta[1:(q+1)] + b
  }

  return(vvec[(p+1):(Tob-q+p)])
}

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
