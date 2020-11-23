#' Estimating an ARCH(q) model
#'
#' @param r data input
#' @param q ARCH order
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#' @export

  arch <- function(r, q, max_iter = 30000, crit = 1e-12) {
    r <- as.matrix(r)

    Tob <- nrow(r)

    r2 <- r^2 # squared residuals
    epsilon <- r2[-c(1:q),] # Arch process

    # generating intital parameters
    theta <- as.matrix(rep(1, q + 1) * 0.5/q)
    m_r2 <- mean(r2)
    theta[1,1] <- m_r2*0.5 # inital value for constant

    Z = YLagCr(r2, q) # Generate regressor matrix

    parameter <- BHHH_arch(r2, q, theta, epsilon, Z, Tob, max_iter = max_iter, crit)

    theta <- parameter[1:(q + 1)]

    scores <- score(epsilon, Z, theta)
    cov_mat <- solve(tcrossprod(scores) / (Tob - q)) / (Tob - q)

    se <- sqrt(diag(cov_mat))

    t_values <- theta/se

    # calculate implied standard deviations
    Z_pre_sample <- YLagCr0(r2, Tob, q , m_r2)

    imp_var_pre_sample <- Z_pre_sample %*% theta
    imp_var <- Z %*% theta
    imp_var <- sqrt(rbind(imp_var_pre_sample, imp_var))

    resid <- r/imp_var

    return(list(theta = theta,
                loglik = parameter[length(parameter)] *(-1),
                se = se,
                t_values = t_values,
                residuals = resid))
}


