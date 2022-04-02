#' Backtesting via Value-at-Risk (VaR)
#'
#' @description Method for calculating VaR from estimated covariance processes (\link{bekk_fit}) or predicted covariances (\link{bekk_forecast}).
#'
#' @param x An object of class "bekkFit" from the function \link{bekk_fit} or an object of class "bekkSpec".
#' @param data A time series if x is of class "bekkSpec".
#' @param window_length An integer specifying the length of the rolling window.
#' @param p A numerical value that determines the confidence level. The default value is set at 0.99 in accordance with the Basel Regulation.
#' @param portfolio_weights A vector determining the portfolio weights to calculate the portfolio VaR. If set to "NULL", the univariate VaR for each series are calculated.
#' @return  Returns a S3 class "var" object containing the VaR forecast and respective confidence bands.
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' # single VaRs of series
#' x2 <- VaR(x1)
#' plot(x2)
#'
#' # VaR of equally-weighted portfolio
#' portfolio_weights <- c(0.5, 0.5)
#' x3 <- VaR(x1, portfolio_weights = portfolio_weights)
#' plot(x3)
#'
#' # VaR of traditional 30/70 weighted bond and stock portfolio
#' portfolio_weights <- c(0.3, 0.7)
#' x4 <- VaR(x1, portfolio_weights = portfolio_weights)
#' plot(x4)
#'
#' }
#'
#' @import xts
#' @import stats
#' @import GAS
#' @export

backtest<- function(x, data=NULL, window_length = 500, p = 0.95, portfolio_weights = NULL,  n.ahead = 1) {
  UseMethod('backtest')
}

#' @export
backtest.bekkFit <-  function(x, data=NULL, window_length = 500, p = 0.95, portfolio_weights = NULL, n.ahead = 1)
{
  data <- x$data
  n <- nrow(data)
  N <- ncol(data)
  n_out = n - window_length





  #portfolio_weights = matrix(portfolio_weights, ncol = N, nrow = 1)
  #out_sample_returns <-  x$data[(window_length+1):n,] %*% t(portfolio_weights)
  if(window_length < 500){
    stop("The supplied window_length must be larger than 500.")
  }
 if(window_length >= n){
   stop("The supplied window_length exeeds the length of the data.")
 }

  if (is.null(portfolio_weights)) {
    hit_rate = numeric(N)
    out_sample_returns <-  x$data[(window_length+1):n,]

    VaR <- matrix(NA, nrow = n_out, ncol = N)

    for(i in 1:n_out){
      spec = bekk_spec()
      fit <- bekk_fit(spec, data[i:(window_length-1+i),])
      forecast <- bekk_forecast(fit, n.ahead = n.ahead, ci = 0.5)
      VaR[i,] = as.matrix(VaR(forecast, p = p, portfolio_weights = portfolio_weights)$VaR[(window_length+1),])
      for(j in 1:N){
        if(VaR[i,j] > out_sample_returns[i,j]){
          hit_rate[j] = hit_rate[j]  + 1
      }

      }

    }
    hit_rate = hit_rate/n_out

    for (i in 1:N) {
      colnames(VaR)[i] <- paste('VaR of', colnames(x$data)[i])
    }
  } else {
    hit_rate = 0
    VaR <- matrix(NA, nrow = nrow(x$data), ncol = 1)
    for(i in 1:n_out){
      spec = bekk_spec()
      fit <- bekk_fit(spec, data[i:(window_length-1+i),])
      forecast <- bekk_forecast(fit, n.ahead = n.ahead, ci = 0.5)
      VaR[i,] = as.matrix(VaR(forecast, p = p, portfolio_weights = portfolio_weights)$VaR[(window_length+1),])
      for(j in 1:N){
        if(VaR[i,j] > out_sample_returns[i,j]){
          hit_rate[j] = hit_rate[j]  + 1
        }

      }

    }
    hit_rate = hit_rate/n_out
    VaR <- as.data.frame(VaR)
  }

  if (inherits(x$data, "ts")) {
    VaR <- ts(VaR, start = time(x$data)[1], frequency = frequency(x$data))
  }

  backtests = suppressWarnings(BacktestVaR(out_sample_returns, VaR, alpha = 1- p))

  result=list(
    VaR = VaR,
    out_sample_returns = out_sample_returns,
    hit_rate = hit_rate,
    backtests = backtests,
    portfolio_weights = portfolio_weights
  )
  class(result) <- c('backtest')
  return(result)
}

