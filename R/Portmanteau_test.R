
#' Performing a Portmanteau test checking for remaining correlation in the empirical co-variances of the estimated BEKK residuals.
#'
#' @description Method for a Portmanteau test of the null hypothesis of no remaining correlation in the co-variances of the estimated BEKK residuals.
#'
#' @param x An object of class "bekkFit" from function \link{bekk_fit}.
#' @param lags Either an integer vector or scalar defining the lag length.
#' @return  Returns a matrix containing the p-values and test statistics.
#'
#' @details Here, the multivariate Portmanteau test of Hosking (1980) is implemented.
#'
#' @references  J. R. M. Hosking (1980). The Multivariate Portmanteau Statistic, Journal of the American Statistical Association, 75:371, 602-608.


#' @import xts
#' @import stats
#' @import ks
#' @export
portmanteau.test <- function(x, lags = 5){
  if(!is.numeric(lags)){
    stop("Please provide a numeric object or vector for 'lags'.")
  }
  if(any(lags<3)){
    stop("Please provide 'lags' larger than 2.")
  }
  UseMethod("portmanteau.test")
}

#' @export
portmanteau.test.bekkFit <- function(x, lags = 5){
  e <- x$e_t
  n <- nrow(e)
  N <- ncol(e)
  #e <- matrix(e, nrow = n, ncol = N)
  e2 <- matrix(NA,nrow = n, ncol = N*(N+1)/2)

  for(i in 1:n){
    e2[i,] <- ks::vech(crossprod(t(e[i,]),e[i,]))
  }
  e=e2

  c_hat <- function(j){
        c= t(e[(j+1):n,]) %*% e[1:(n-j),]
        return(c/n)
  }
  c_0 = c_hat(0)
  c_0_inv = solve(c_0)

  Q <- function(lgs){
    q=0
   for(i in 1:lgs){
     c_temp = t(c_hat(i))
     q=q+sum(diag(c_temp%*%c_0_inv%*%c_temp%*%c_0_inv))
   }
    return(q)
  }
  p_val_q <- function(c, lgs){
    return(1-pchisq(c, df=(lgs-2)*(N^2)))
  }

  summary_res <- function(lgs){
    t_stat = Q(lgs)
    return(c(lgs, t_stat, (lgs-2)*(N^2), p_val_q(t_stat,lgs)))
  }


  res = t(sapply(lags,FUN=summary_res))
  dimnames(res) <- list(rep("", length(lags)),c("lags","statistic","df","p-value"))
  cat(paste("\n", "Portmanteau Test", "\n", sep = ""))
  underScore <- paste(rep("-", nchar("Portmanteau Test")), collapse = "")
  cat(underScore, "\n")
  return(res)
}
