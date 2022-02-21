#' Estimating multivariate volatility impulse response functions (VIRF) for BEKK models
#'
#' @description Method for estimating VIRFs of N-dimensional BEKK models.
#'
#' @param fit An object of class "bekkfit" from function \link{bekk_fit}.
#' @param time Time instace to calculate VIRFs for.
#' @param q Quantiles for VIRFs.
#' @param periods Periods ahead for VIRFs.
#' @return  Returns an object of class "bekkVIRF".
#'
#' @import xts
#' @import stats
#' @export

bekk_virf <- function(fit ,time = 1, q, periods = 10){

  if (!inherits(fit, 'bekkfit')) {
    stop('Please provide and object of class "bekkFit" for fit')
  }

  UseMethod('bekk_virf')

}

#' @export
bekk_virf.bekk <- function(fit, time, q) {


  N <- ncol(fit$data)
  q_len <- length(q)
  H <- fit$H_t[[time]]
  #get quantiles of returns
  residuals = fit$resiudals

  quantiles = quantile(residuals[,1],probs=q[1])
  VIRF =

  result <- list()
  class(result) <- c('bekkFit', 'bekk')
  return(result)
}
