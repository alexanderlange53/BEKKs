#' Estimating multivariate volatility impulse response functions (VIRF) for BEKK models
#'
#' @description Method for estimating VIRFs of N-dimensional BEKK models.
#'
#' @param fit An object of class "bekkfit" from function \link{bekk_fit}.
#' @param time Time instace to calculate VIRFs for.
#' @param q A vector specifying the quantiles to be considered for a shock on which basis the VIRFs are generated.
#' @param periods An integer defining the number periods for which the VIRFs are generated.
#' @return  Returns an object of class "bekkVIRF".
#'
#' @import xts
#' @import stats
#' @export

virf <- function(fit ,time = 1, q = 0.05, index_series=1, periods = 10){

  if (!inherits(fit, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for fit')
  }

  UseMethod('bekk_virf')

}

#' @export
virf.bekk <- function(fit, time = 1, q = 0.05, index_series=1, periods = 10) {

  N <- ncol(fit$data)
  H <- matrix(fit$H_t[time,],N,N)
  #get quantiles of returns
  residuals = fit$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))


  VIRF = virf_bekk(H, fit$A, fit$G, matrix(shocks,ncol=N, nrow = 1), periods)

  result <- list(VIRF=VIRF,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 fit=fit)
  class(result) <- c('bekkVirf','bekkFit', 'bekk')
  return(result)
}

virf.bekka <- function(fit, time = 1, q = 0.05, index_series=1, periods = 10) {

  N <- ncol(fit$data)
  H <- matrix(fit$H_t[time,],N,N)
  e
  #get quantiles of returns
  residuals = fit$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))


  VIRF =  virf_bekka(H_t, fit$A, fit$B, fit$G, fit$signs, fit$expected_signs, shocks, periods)


  result <- list(VIRF=VIRF,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 fit=fit)
  class(result) <- c('bekkVirf','bekkFit', 'bekka')
  return(result)
}
