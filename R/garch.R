#' Estimating an GARCH(p, q) model
#'
#' @param r data input
#' @param p GARCH order
#' @param q ARCH order
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#' @export

  garch <- function(r, p, q, max_iter = 300000, crit = 1e-12) {
    r <- as.matrix(r)

    Tob <- nrow(r)

    r2 <- r^2 # squared residuals
    epsilon2 <- r2[-c(1:q),] # GARCH process
    ucvar <- sum(r2) / (Tob - q) # unconditional variance

    # generating intital parameters
    theta <- as.matrix(c(ucvar*0.05, rep(1, q) * 0.05/q, rep(1, p) * 0.9/p))

    Z = YLagCr(r2, q) # Generate regressor matrix

    parameter <- BHHH_garch(r2, q, p, theta, epsilon2, Z, Tob, max_iter, crit, ucvar)

    theta <- parameter[1:(q + p + 1)]

    scores <- score_garch(epsilon2, Z, Tob, q, p, theta, ucvar)
    cov_mat <- solve(tcrossprod(scores) / (Tob - q)) / (Tob - q)

    se <- sqrt(diag(cov_mat))

    t_values <- theta/se

    # calculate implied standard deviations
    Z_pre_sample <- YLagCr0(r2, Tob, q , ucvar)
    Z_pre_sample <- rbind(Z_pre_sample, Z)

    imp_var <- greatZ0_garch(Tob, Z_pre_sample, q, p, ucvar, theta)

    resid <- r/imp_var

    return(list(theta = theta,
                loglik = parameter[length(parameter)] *(-1),
                se = se,
                t_values = t_values,
                residuals = resid))
}


