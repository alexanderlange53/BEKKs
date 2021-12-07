#' Forecasting conditional volatilities with BEKK models
#'
#' @param x A fitted bekk model of class bekk from the \link{bekk} function
#' @param n.ahead Number of periods to forecast conditional volatility. Default is a one-period ahead forecast.
#'
#' @examples
#' \donttest{
#'
#' data(StocskBonds)
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)'
#'
#' x2 <- bekk_forecast(x1, n.ahead = 1)
#'
#' }
#' @export
bekk_forecast <- function(x, n.ahead = 1){

  if (!inherits(x, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for "x".')
  }


  UseMethod('bekk_forecast')

}

#' @export
bekk_forecast.bekk <- function(x, n.ahead = 1) {
  N <- ncol(x$data)
  NoBs <- nrow(x$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t <- vector(mode = "list",length = n.ahead+1)
  H_t[[1]] <- matrix(x$H_t[NoBs,],nrow = N, ncol = N)
  current_returns <- t(x$data[NoBs,])

  for(i in 1:n.ahead){
    H_t[[i+1]] <- t(x$C0) %*% x$C0 + t(x$A) %*% t(current_returns) %*% current_returns %*% x$A + t(x$G) %*% H_t[[i]] %*% x$G
    current_returns <- eigen_value_decomposition(H_t[[i+1]])
  }

  sigma_t <- matrix(NA, nrow = n.ahead, ncol = N^2)
  for (i in 2:(n.ahead+1)){
      tm2 <- sqrt(solve(diag(diag(H_t[[i]]))))%*%H_t[[i]]%*%sqrt(solve(diag(diag(H_t[[i]]))))
      diag(tm2) <- sqrt(diag(H_t[[i]]))
      sigma_t[i-1,] <- c(tm2)
  }

  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  H_t_f <- matrix(NA, nrow = n.ahead, ncol = ncol(x$data)^2)

  for (i in 2:(n.ahead+1)){
    H_t_f[i-1, ] <- c(H_t[[i]])
  }

  result <- list(
    volatility_forecast = sigma_t,
    H_t_forecast = H_t_f,
    n.ahead = n.ahead,
    bekkfit = x
  )
  class(result) <- c('bekkForecast', 'bekk')
  return(result)
}

#' @export
bekk_forecast.bekka <- function(x, n.ahead = 1) {
  N <- ncol(x$data)
  NoBs <- nrow(x$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t = vector(mode = "list",length=n.ahead+1)
  H_t[[1]] = matrix(x$H_t[NoBs,],nrow = N, ncol = N)
  current_returns = t(x$data[NoBs,])
  H_t[[2]]=t(x$C0) %*% x$C0 + t(x$A) %*% t(current_returns) %*% current_returns %*% x$A + indicatorFunction(as.matrix(current_returns), x$signs) * t(x$B) %*% t(current_returns) %*% current_returns %*% x$B + t(x$G) %*% H_t[[1]] %*% x$G
  #Some work left. Expected signs need to be changed 
  expected_signs=0
  for(i in 2:n.ahead){
    H_t[[i+1]]=t(x$C0) %*% x$C0 + t(x$A) %*% H_t[[i+1]]s %*% x$A + expected_signs * t(x$B) %*% H_t[[i+1]] %*% x$B + t(x$G) %*% H_t[[i]] %*% x$G
  
  }

  sigma_t = matrix(NA, nrow = n.ahead, ncol = N^2)
  for (i in 1:(n.ahead+1)){
    tm2 <- sqrt(solve(diag(diag(H_t[[i]]))))%*%H_t[[i]]%*%sqrt(solve(diag(diag(H_t[[i]]))))
    diag(tm2) <- sqrt(diag(H_t[[i]]))
    sigma_t[i-1,] <- c(tm2)
  }

  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]


  H_t_f <- matrix(NA, nrow = n.ahead, ncol = ncol(x$data)^2)

  for (i in 1:(n.ahead+1)){
    H_t_f[i-1, ] <- c(H_t[[i]])
  }

  result <- list(
    volatility_forecast = sigma_t,
    H_t_forecast = H_t_f,
    n.ahead = n.ahead,
    bekkfit = x
  )
  class(result) <-  c('bekkForecast', 'bekka')
  return(result)
}
