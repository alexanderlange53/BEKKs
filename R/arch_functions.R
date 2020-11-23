# Script containing the functions to estimate the ARCH(q) parameter

BHHH_garch <- function(r2, q, theta, epsilon2, Z, Tob, max_iter, crit){
  # BHHH algorithm based on first order derivatives and outer product for ARCH(q) model

  step_width <- seq(1, 0, length.out = 100) # find the optimal step width

  quit <- 0
  iter <- 1

  while (quit == 0 & iter < max_iter) {

    score_function <- score(epsilon2, Z, theta)

    out_prod_grad <- score_function %*% t(score_function)

    sum_grad <- rowSums(score_function)

    theta_hat <- matrix(0, nrow = nrow(theta), ncol = 100)
    for (i in 1:100) {
      theta_hat[,i] <- theta + step_width[i] * solve(out_prod_grad) %*% sum_grad
    }

    theta_hat <- abs(theta_hat)

    candidates <- matrix(0, 100, 1)

    candidates[1, ] <- likelihood_arch(Z, Tob, q, theta, epsilon2)

    quit_inner <- 0
    iter_inner <- 2

    while (iter_inner < 100 & quit_inner == 0) {
      candidates[iter_inner, ] <- likelihood_arch(Z, Tob, q, theta_hat[, iter_inner], epsilon2)

      if (candidates[iter_inner] < candidates[iter_inner-1]) {
        quit_inner <- 1
      }
      iter_inner <- iter_inner + 1
    }

    min_lik <- which.min(candidates)

    best_candidate <- candidates[min_lik,]

    if ((best_candidate - candidates[1,])^2 / abs(candidates[1,]) <  crit) {
      quit <- 1
    }

      theta <- as.matrix(theta_hat[, min_lik])
    iter <- iter + 1
  }

  lik_optim <- likelihood_arch(Z, Tob, q, theta, epsilon2)

  return(c(theta, lik_optim))
}

score <- function (epsilon2, Z, theta) {
  # Calculates the first order derivatives

  sigma2 <- Z %*% theta

  sigg1 <- kronecker(sigma2^(-1), matrix(1, nrow = 1, ncol(Z))) * Z
  sigg2 <- kronecker(sigma2^(-2), matrix(1, nrow = 1, ncol(Z))) * Z

  lt <- 0.5 * ((-1) * sigg1 + epsilon2 * sigg2)

  return(t(lt))
}

score_garch <- function(epsilon2, Z, Tob, q, p, theta, ucvar) {
  vvec <- greatZ_garch(Tob, Z, q, p, ucvar, theta)
  vvec2 <- c(ucvar, ucvar, vvec)
  vvec2 <- vvec2[1:(Tob-q)]
  vvec1 <- c(ucvar, vvec)
  vvec1 <- vvec1[1:(Tob-q)]

  aa <- ucvar/(1 - theta[(q + 2):nrow(theta)] %*% rep(1, nrow(theta) - q -q))
  aa2 <- 1/(1 - theta[(q+2):nrow(theta)] %*% rep(1, nrow(theta) - q -q))

  gz <- cbind(Z, vvec1, vvec2)
  g <- nrow(theta)

  gz <- gz[,1:g]
  vabl <- matrix(0, nrow=(Tob-q), ncol = g)

  e <- aa %*% matrix(1, p, q)
  vabl <- c(e, vabl)
}

likelihood_arch <- function(Z, Tob, q, p, theta, epsilon2) {
  # likelihood function for specific parameter 'theta'

  #sigma2 <- Z %*% theta
  sigma2 <- greatZ_garch(Tob, Z, q, p, ucvar, theta)

  term2 <- epsilon2 / sigma2

  loglik <- 0.5 * (Tob - q) * log(2 * pi) + 0.5 * sum(log(sigma2)) + 0.5 * sum(term2)
  return(loglik)
}

# YLagCr0A <- function(r2, Tob, q, m_r2) {
#
#   help_m <- matrix(1, nrow = nrow(r2) + q, ncol = q +1)
#   for (j in 2:(q + 1)) {
#     k <- 1
#     for (i in 1:q) {
#       if (i < j) {
#         help_m[i, j] <- m_r2
#       } else {
#         help_m[i, j] <- r2[k,]
#         k <- k + 1
#       }
#     }
#   }
#
#   return(help_m[1:q, 1:(q+1)])
# }

# Generates vector of variances
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
