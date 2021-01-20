
library(maxLik)
bigarchR = function(theta, returns) {
    # determine C
    A      = matrix(theta[1:4], 2)
    G      = matrix(theta[5:8], 2)
    sigma0 = (t(returns) %*% returns) / nrow(returns)
    C2v    = (diag(4) - t(kronecker(A, A)) - t(kronecker(G, G))) %*% c(sigma0)
    C      = chol(matrix(C2v, 2))
    
    val    = maxLik(
        logLik  = loglike_bekk,
        start   = c(c(C)[-2], theta),
        method  = "BHHH",
        control = list(
            iterlim    = 1e9,
            printLevel = 0,
            steptol    = 1e-12,
            tol        = 1e-10
        ),
        data    = returns
    )
    
    # return estimated coefficients and loglikelihood value
    return(list(coeff  = val$estimate,
                maxlik = val$maximum))
}


# validate coefficients (constraints)
valid_bekk <- function(C, A, G) {
    # condition for positive-definit covariance
    if (any(diag(C) < 0)) {
        return(FALSE)
    }
    
    # condition for uniqueness
    if (A[1, 1] < 0 || G[1, 1] < 0) {
        return(FALSE)
    }
    
    # check stationarity for BEKK(1,1): Engle & Kroner (1995), Prop. 2.7
    if (!all(abs(eigen(kronecker(A, A)
                       + kronecker(G, G))$values) < 1)) {
        return(FALSE)
    }
    
    return(TRUE)
}


# compute (co)variance processes
comph_bekk <- function(C, A, G, data) {
    # dimensions
    n      = nrow(data)
    
    # redefine for convenience
    CC     = t(C) %*% C
    At     = t(A)
    Gt     = t(G)
    
    # compute uncond. covariance matrix
    Hu     = (t(data) %*% data) / nrow(data)
    
    # compute H trajectory
    H      = vector(mode = "list", n)
    H[[1]] = Hu
    for (i in 2:n) {
        H[[i]] <- (CC + At %*% (data[i - 1,] %*% t(data[i - 1,])) %*% A
                   + Gt %*% H[[i - 1]] %*% G)
    }
    return(H)
}


# compute loglikelihood vector
loglike_bekk <- function(theta, data) {
    # convert to matrices
    C = matrix(c(theta[1], 0, theta[2:3]), 2)
    A = matrix(theta[4:7], 2)
    G = matrix(theta[8:11], 2)
    
    # check constraints
    if (!valid_bekk(C, A, G)) {
        return(NA)
    }
    
    # compute H
    H   = comph_bekk(C, A, G, data)
    
    # compute llv
    n   = nrow(data)
    llv = numeric(n)
    for (i in 1:n) {
        llv[i] = c(t(data[i,]) %*% solve(H[[i]]) %*% data[i,])
    }
    llv = -0.5 * (2 * log(2 * pi) + log(sapply(H, det)) + llv)
    return(llv)
}
