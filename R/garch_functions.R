# Script containing the functions to estimate the GARCH(p, q) parameter

BHHH_garch <- function(r2, q, p, theta, epsilon2, Z, Tob, max_iter, crit, ucvar){
  # BHHH algorithm based on first order derivatives and outer product for ARCH(q) model

  step_width <- seq(1, 0, length.out = 100) # find the optimal step width

  quit <- 0
  iter <- 1

  while (quit == 0 & iter < max_iter) {

    theta <- matrix(theta, nrow = p+q+1, ncol = 1)
    score_function <- ScoreGarch(matrix(epsilon2, 1, length(epsilon2)), Z, Tob, q, p, theta, ucvar)

    out_prod_grad <- score_function %*% t(score_function)

    sum_grad <- rowSums(score_function)

    theta_hat <- matrix(0, nrow = nrow(theta), ncol = 100)
    for (i in 1:100) {
      theta_hat[,i] <- theta + step_width[i] * solve(out_prod_grad) %*% sum_grad
    }

    theta_hat <- abs(theta_hat)

    candidates <- matrix(0, 100, 1)

    candidates[1, ] <- LikelihoodGarch(Z, Tob, q, p, matrix(theta, ncol = 1), matrix(epsilon2, 1, length(epsilon2)), ucvar)

    quit_inner <- 0
    iter_inner <- 2

    while (iter_inner < 100 & quit_inner == 0) {
      candidates[iter_inner, ] <- LikelihoodGarch(Z, Tob, q, p, matrix(theta_hat[, iter_inner], ncol= 1), matrix(epsilon2, 1, length(epsilon2)), ucvar)
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

  lik_optim <- LikelihoodGarch(Z, Tob, q, p, matrix(theta, ncol = 1), matrix(epsilon2, 1, length(epsilon2)), ucvar)

  return(c(theta, lik_optim))
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
