#' Forecasting conditional volatilities with BEKK models
#'
#' @param spec A fitted bekk model of class bekk from the \link{bekk} function
#' @param forecast_length Number of periods to forecast conditional volatility. Default is a one-period ahead forecast.
#'
#' @examples
#' \donttest{
#'
#' data(bivariate)
#' x1 <- bekk_fit(BI, init_values = NULL,
#' QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' }
#' @export
bekk_forecast <- function(spec, forecast_length = 1){

  if (!inherits(spec, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for spec.')
  }


  UseMethod('bekk_forecast')

}

#' @export
bekk_forecast.bekk <- function(spec, forecast_length = 1) {
  xx <- spec
  n <- ncol(xx$data)
  NoBs <- nrow(xx$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t = vector(mode = "list",length=forecast_length+1)
  H_t[[1]] = matrix(xx$H_t[NoBs,],nrow = n, ncol = n)
  current_returns = xx$data[NoBs,]

  for(i in 1:forecast_length){
    H_t[[i+1]]=t(xx$C0) %*% xx$C0 + t(xx$A) %*% t(current_returns) %*% current_returns %*% xx$A + t(xx$G) %*% H_t[[i]] %*% xx$G
    current_returns = t(as.matrix(rnorm(n))) %*% eigen_value_decomposition(H_t[[i+1]])
  }

  sigma_t = matrix(NA, nrow = forecast_length, ncol = n^2)
  for (i in 1: forecast_length){
      tm2 <- sqrt(solve(diag(diag(H_t[[i]]))))%*%H_t[[i]]%*%sqrt(solve(diag(diag(H_t[[i]]))))
      diag(tm2) <- sqrt(diag(H_t[[i]]))
      sigma_t[i,] <- c(tm2)
  }
  #Hier in Zukunft noch VaR Forecasts?



  result <- list(
    volatility_forecast = sigma_t
  )
  class(result) <- c('bekkFit', 'bekk', 'forecast')
  return(result)
}

#' @export
bekk_forecast.bekka <- function(spec, forecast_length = 1) {
  xx <- spec
  n <- ncol(xx$data)
  NoBs <- nrow(xx$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t = vector(mode = "list",length=forecast_length+1)
  H_t[[1]] = matrix(xx$H_t[NoBs,],nrow = n, ncol = n)
  current_returns = xx$data[NoBs,]

  for(i in 1:forecast_length){
    H_t[[i+1]]=t(xx$C0) %*% xx$C0 + t(xx$A) %*% t(current_returns) %*% current_returns %*% xx$A + indicatorFunction(as.matrix(current_returns),xx$signs) * t(xx$B) %*% t(current_returns) %*% current_returns %*% xx$B + t(xx$G) %*% H_t[[i]] %*% xx$G
    current_returns = t(as.matrix(rnorm(n))) %*% eigen_value_decomposition(H_t[[i+1]])
  }

  sigma_t = matrix(NA, nrow = forecast_length, ncol = n^2)
  for (i in 1: forecast_length){
    tm2 <- sqrt(solve(diag(diag(H_t[[i]]))))%*%H_t[[i]]%*%sqrt(solve(diag(diag(H_t[[i]]))))
    diag(tm2) <- sqrt(diag(H_t[[i]]))
    sigma_t[i,] <- c(tm2)
  }
  #Hier in Zukunft noch VaR Forecasts?



  result <- list(
    volatility_forecast = sigma_t
  )
  class(result) <- c('bekkFit', 'bekk', 'forecast')
  return(result)
}
