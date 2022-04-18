#' Backtesting via Value-at-Risk (VaR)
#'
#' @description Method for calculating VaR from estimated covariance processes (\link{bekk_fit}).
#'
#' @param x An object of class "bekkFit" from the function \link{bekk_fit}.
#' @param window_length An integer specifying the length of the rolling window.
#' @param p A numerical value that determines the confidence level. The default value is set at 0.99 in accordance with the Basel Regulation.
#' @param portfolio_weights A vector determining the portfolio weights to calculate the portfolio VaR. If set to "NULL", the univariate VaR for each series are calculated.
#' @param n.ahead Number of periods to forecast conditional volatility. Default is a one-period ahead forecast.
#' @return  Returns a S3 class "backtest" object containing the VaR forecast, out-of-sample returns and backtest statistics according to the R-package "GAS". conf
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' # backtesting
#' x2 <- backtest(x1, window_length = 6000)
#' plot(x2)
#'
#' }
#'
#' @import xts
#' @import stats
#' @importFrom GAS BacktestVaR
#' @importFrom lubridate date_decimal
#' @export

backtest<- function(x, window_length = 500, p = 0.95, portfolio_weights = NULL,  n.ahead = 1) {
  UseMethod('backtest')
}

#' @export
backtest.bekkFit <-  function(x, window_length = 500, p = 0.95, portfolio_weights = NULL, n.ahead = 1)
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

  if((n-window_length) < n.ahead){
    stop("The supplied 'n.ahead' exceeds the forecasting horizon.")
  }
  if (is.null(portfolio_weights)) {
    hit_rate = numeric(N)
    out_sample_returns <-  x$data[(window_length+1):n,]

    VaR <- matrix(NA, nrow = n_out, ncol = N)

    i = 1
    while(i <= n_out){
      spec = bekk_spec()
      fit <- bekk_fit(spec, data[i:(window_length-1+i),])
      forecast <- bekk_forecast(fit, n.ahead = n.ahead, ci = 0.5)
      VaR[i:(i+n.ahead-1),] = as.matrix(VaR(forecast, p = p, portfolio_weights = portfolio_weights)$VaR[(window_length+1):(window_length+n.ahead),])


      for(j in 1:N){
      hit_rate[j]= hit_rate[j]  + sum(VaR[i:(i+n.ahead-1),j] > out_sample_returns[i:(i+n.ahead-1),j])
      }

      if(n.ahead > 1 && i >= (n_out-n.ahead)){
        n.ahead = 1
      }
      i = i + n.ahead

    }
    hit_rate = hit_rate/n_out
    backtests = list()



    VaR <- as.data.frame(VaR)
    for (i in 1:N) {
      backtests[[i]] = suppressWarnings(BacktestVaR(out_sample_returns[,i], VaR[,i], alpha = 1- p))
      colnames(VaR)[i] <- paste('VaR of', colnames(x$data)[i])
    }
  } else {
    out_sample_returns = x$data[(window_length+1):n,] %*% matrix(portfolio_weights, ncol = 1, nrow = N)
    hit_rate = 0

    VaR <- matrix(NA, nrow = n_out, ncol = 1)
    i = 1
    while(i <= n_out){

      spec = x$spec
      fit <- bekk_fit(spec, data[i:(window_length-1+i),])
      forecast <- bekk_forecast(fit, n.ahead = n.ahead, ci = 0.5)
      VaR[i:(i+n.ahead-1),] = as.matrix(VaR(forecast, p = p, portfolio_weights = portfolio_weights)$VaR[(window_length+1):(window_length+n.ahead),])


      hit_rate= hit_rate  + sum(VaR[i:(i+n.ahead-1),] > out_sample_returns[i:(i+n.ahead-1),])


     if(n.ahead > 1 && i >= (n_out-n.ahead)){
        n.ahead = 1
     }
       i = i + n.ahead

    }
    hit_rate = hit_rate/n_out
    backtests= suppressWarnings(GAS::BacktestVaR(out_sample_returns, VaR, alpha = 1- p))
    VaR <- as.data.frame(VaR)
  }
  out_sample_returns = as.data.frame(out_sample_returns)

  if (inherits(x$data, "ts")) {
    VaR <- xts(VaR, order.by = date_decimal(time(x$data[(window_length+1):n,])))
    out_sample_returns <- xts(out_sample_returns, order.by = date_decimal(time(x$data[(window_length+1):n,])))
  }else if(inherits(x$data, "xts") || inherits(x$data, "zoo") ){
    VaR <- xts(VaR, order.by = time(x$data[(window_length+1):n,]))
    out_sample_returns <- xts(out_sample_returns, order.by = time(x$data[(window_length+1):n,]))

  }



  result=list(
    VaR = VaR,
    out_sample_returns = out_sample_returns,
    hit_rate = hit_rate,
    backtests = backtests,
    portfolio_weights = portfolio_weights
  )
  class(result) <- c('backtest', 'bekkFit')
  return(result)
}

