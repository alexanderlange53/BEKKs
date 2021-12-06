#' Calculating Value-at-Risk (VaR)
#'
#' @description Method for calculating VaR from estimated covariance processes (\link{bekk_fit}) or predicted covariances (\link{bekk_forecast}).
#'
#' @param x An object of class "bekkFit" from the function \link{bekk_fit} or an object of class "bekkForecast" from the function \link{bekk_forecast}.
#' @param p A numerical value that determines the confidence level. The default value is set at 0.99 in accordance with the Basel Regulation.
#' @param portfolio_weights A vector determing the portfolio weights to calculate the portfolio VaR. If set to "NULL", the univariate VaR for each series are calculated.
#'
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
#' @export

VaR <- function(x, p = 0.99, portfolio_weights = NULL) {
  UseMethod('VaR')
}

#' @export
VaR.bekkFit <-  function(x, p = 0.99, portfolio_weights = NULL)
{
  alpha = p

  if (is.null(portfolio_weights)) {
    columns = ncol(x$data)
    csd <- extract_csd(x)
    VaR <- matrix(NA, nrow = nrow(x$data), ncol = ncol(x$data))

    for(column in 1:columns) {
      r = as.vector(na.omit(x$data[,column]))
      if (!is.numeric(r)) stop("The selected column is not numeric")
      m2 =  csd[, column]
      VaR[, column] = - mean(r) - qnorm(alpha)*m2
      VaR <- as.data.frame(VaR)

      for (i in 1:ncol(x$data)) {
        colnames(VaR)[i] <- paste('VaR of', colnames(x$data)[i])
      }
    }
  } else {
    VaR <- matrix(NA, nrow = nrow(x$data), ncol = 1)
    for(i in 1:nrow(x$H_t)) {
      VaR[i,] <- -qnorm(alpha)*portfolio_weights%*%eigen_value_decomposition(matrix(x$H_t[i,], ncol = ncol(x$data)))%*%portfolio_weights
    }
    VaR <- as.data.frame(VaR)
  }

  if (inherits(x$data, "ts")) {
    VaR <- ts(VaR, start = time(x$data)[1], frequency = frequency(x$data))
  }

  result <- list(VaR = VaR,
                 p = p,
                 portfolio_weights = portfolio_weights,
                 bekk = x)
  class(result) <- 'var'
  return(result)
}
