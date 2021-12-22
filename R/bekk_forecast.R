#' Forecasting conditional volatilities with BEKK models
#'
#' @param x A fitted bekk model of class bekk from the \link{bekk_fit} function
#' @param n.ahead Number of periods to forecast conditional volatility. Default is a one-period ahead forecast.
#' @param ci Floating point in [0,1] defining the niveau for confidence bands of the conditional volatility forecast. Default is 95% niveau confidence bands.
#'
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)'
#'
#' x2 <- bekk_forecast(x1, n.ahead = 1)
#'
#' }
#' @export
bekk_forecast <- function(x, n.ahead = 1, ci = 0.95){

  if (!inherits(x, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for "x".')
  }


  UseMethod('bekk_forecast')

}

#' @export
bekk_forecast.bekk <- function(x, n.ahead = 1, ci = 0.95) {
  N <- ncol(x$data)
  NoBs <- nrow(x$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t <- vector(mode = "list",length = n.ahead+1)
  H_t[[1]] <- matrix(x$H_t[NoBs,],nrow = N, ncol = N)
  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)

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

  # Generating confidence intervals

  score_final = score_bekk(x$theta, x$data)
  s1_temp = diag(solve(t(score_final) %*% score_final),names=T)
  s1 = sqrt(s1_temp)

  lower_theta = x$theta - qnorm(ci)*s1
  upper_theta = x$theta + qnorm(ci)*s1

  H_t_lower <- vector(mode = "list",length = n.ahead+1)
  H_t_upper <- vector(mode = "list",length = n.ahead+1)
  x_lower <- coef_mat(lower_theta, N)
  x_upper <- coef_mat(upper_theta, N)
  #check if t(x_lower$C0) or x_lower$C0 for sigma_bekk
  H_t_lower[[1]] <- matrix(sigma_bekk(x$data, t(x_lower$c0), x_lower$a, x_lower$g)$sigma_t[NoBs,], nrow = N, ncol = N)
  H_t_upper[[1]] <- matrix(sigma_bekk(x$data,t(x_upper$c0),x_upper$a, x_upper$g)$sigma_t[NoBs,], nrow = N, ncol = N)
  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)


  for(i in 1:n.ahead){
    H_t_lower[[i+1]] <- t(x_lower$c0) %*% x_lower$c0 + t(x_lower$a) %*% t(current_returns) %*% current_returns %*% x_lower$a + t(x_lower$g) %*% H_t[[i]] %*% x_lower$g
    current_returns <- eigen_value_decomposition(H_t_lower[[i+1]])
  }

  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)

  for(i in 1:n.ahead){
    H_t_upper[[i+1]] <- t(x_upper$c0) %*% x_upper$c0 + t(x_upper$a) %*% t(current_returns) %*% current_returns %*% x_upper$a + t(x_upper$g) %*% H_t[[i]] %*% x_upper$g
    current_returns <- eigen_value_decomposition(H_t_upper[[i+1]])
  }
  sigma_t_lower <- matrix(NA, nrow = n.ahead, ncol = N^2)

  for (i in 2:(n.ahead+1)){
    tm2 <- sqrt(solve(diag(diag(H_t_lower[[i]]))))%*%H_t_lower[[i]]%*%sqrt(solve(diag(diag(H_t_lower[[i]]))))
    diag(tm2) <- sqrt(diag(H_t_lower[[i]]))
    sigma_t_lower[i-1,] <- c(tm2)
  }
  sigma_t_upper <- matrix(NA, nrow = n.ahead, ncol = N^2)
  for (i in 2:(n.ahead+1)){
    tm2 <- sqrt(solve(diag(diag(H_t_upper[[i]]))))%*%H_t_upper[[i]]%*%sqrt(solve(diag(diag(H_t_upper[[i]]))))
    diag(tm2) <- sqrt(diag(H_t_upper[[i]]))
    sigma_t_upper[i-1,] <- c(tm2)
  }


  colnames(sigma_t_lower) <- rep(1, N^2)
  colnames(sigma_t_upper) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t_lower)[k2] <- paste('Lower CI conditional standard deviation of \n', colnames(x$data)[k])
        colnames(sigma_t_upper)[k2] <- paste('Upper CI conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t_lower)[k2] <- paste('Lower CI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        colnames(sigma_t_upper)[k2] <- paste('Upper CI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])

        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t_lower <- sigma_t_lower[, which(colSums(elim) == 1)]
  sigma_t_upper <- sigma_t_upper[, which(colSums(elim) == 1)]

  H_t_f_lower <- H_t_f_upper <- matrix(NA, nrow = n.ahead + 1, ncol = ncol(x$data)^2)

  for (i in 1:(n.ahead+1)){
    H_t_f_lower[i, ] <- c(H_t_lower[[i]])
    H_t_f_upper[i, ] <- c(H_t_upper[[i]])
  }

  if (inherits(x$data, "ts")) {
    sigma_t <- ts(sigma_t, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
    sigma_t_lower <- ts(sigma_t_lower, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
    sigma_t_upper <- ts(sigma_t_upper, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
  }
  else if(inherits(x$data, "xts") || inherits(x$data, "zoo") ){
    sigma_t <- xts(matrix(sigma_t, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                   by = periodicity(x$data)$units))
    sigma_t_lower <- xts(matrix(sigma_t_lower, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                         by = periodicity(x$data)$units))
    sigma_t_upper <- xts(matrix(sigma_t_upper, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                         by = periodicity(x$data)$units))
  }


  result <- list(
    volatility_forecast = sigma_t,
    volatility_lower_conf_band = sigma_t_lower,
    volatility_upper_conf_band = sigma_t_upper,
    H_t_forecast = H_t_f,
    H_t_lower_conf_band = H_t_f_lower,
    H_t_upper_conf_band = H_t_f_upper,
    n.ahead = n.ahead,
    bekkfit = x
  )
  class(result) <- c('bekkForecast', 'bekk')
  return(result)
}

#' @export
bekk_forecast.bekka <- function(x, n.ahead = 1, ci = 0.95) {
  N <- ncol(x$data)
  NoBs <- nrow(x$data)
  #var_process <- sigma_bekk(xx$data, xx$C0, xx$A, xx$G)
  H_t = vector(mode = "list",length=n.ahead+1)
  H_t[[1]] = matrix(x$H_t[NoBs,],nrow = N, ncol = N)
  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)

  H_t[[2]]=t(x$C0) %*% x$C0 + t(x$A) %*% t(current_returns) %*% current_returns %*% x$A + indicatorFunction(as.matrix(current_returns), x$signs) * t(x$B) %*% t(current_returns) %*% current_returns %*% x$B + t(x$G) %*% H_t[[1]] %*% x$G

  expected_signs=expected_indicator_value(x$data,x$signs)
  for(i in 2:n.ahead){
    H_t[[i+1]]=t(x$C0) %*% x$C0 + t(x$A) %*% H_t[[i]] %*% x$A + expected_signs * t(x$B) %*% H_t[[i]] %*% x$B + t(x$G) %*% H_t[[i]] %*% x$G

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


  #95% confidence interval

  score_final = score_asymm_bekk(x$theta, x$data, x$signs)
  s1_temp = diag(solve(t(score_final) %*% score_final),names=T)
  s1 = sqrt(s1_temp)

  lower_theta = x$theta - qnorm(ci)*s1
  upper_theta = x$theta + qnorm(ci)*s1

  H_t_lower <- vector(mode = "list",length = n.ahead+1)
  H_t_upper <- vector(mode = "list",length = n.ahead+1)
  x_lower <- coef_mat_asymm(lower_theta, N)
  x_upper <- coef_mat_asymm(upper_theta, N)
  #check if t(x_lower$C0) or x_lower$C0 for sigma_bekk
  H_t_lower[[1]] <- matrix(sigma_bekk_asymm(x$data, t(x_lower$c0), x_lower$a, x_lower$b, x_lower$g, x$signs)$sigma_t[NoBs,],nrow = N, ncol = N)
  H_t_upper[[1]] <- matrix(sigma_bekk_asymm(x$data, t(x_upper$c0), x_upper$a, x_upper$b, x_upper$g, x$signs)$sigma_t[NoBs,],nrow = N, ncol = N)
  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)


  for(i in 1:n.ahead){
    H_t_lower[[i+1]] <- t(x_lower$c0) %*% x_lower$c0 + t(x_lower$a) %*% t(current_returns) %*% current_returns %*% x_lower$a + expected_signs * t(x_lower$b) %*% H_t_lower[[i]] %*% x_lower$b + t(x_lower$g) %*% H_t[[i]] %*% x_lower$g
    current_returns <- eigen_value_decomposition(H_t_lower[[i+1]])
  }

  current_returns <- matrix(c(x$data[NoBs,]), nrow = 1)

  for(i in 1:n.ahead){
    H_t_upper[[i+1]] <- t(x_upper$c0) %*% x_upper$c0 + t(x_upper$a) %*% t(current_returns) %*% current_returns %*% x_upper$a + expected_signs * t(x_upper$b) %*% H_t_upper[[i]] %*% x_upper$b + t(x_upper$g) %*% H_t[[i]] %*% x_upper$g
    current_returns <- eigen_value_decomposition(H_t_upper[[i+1]])
  }
  sigma_t_lower <- matrix(NA, nrow = n.ahead, ncol = N^2)

  for (i in 2:(n.ahead+1)){
    tm2 <- sqrt(solve(diag(diag(H_t_lower[[i]]))))%*%H_t_lower[[i]]%*%sqrt(solve(diag(diag(H_t_lower[[i]]))))
    diag(tm2) <- sqrt(diag(H_t_lower[[i]]))
    sigma_t_lower[i-1,] <- c(tm2)
  }
  sigma_t_upper <- matrix(NA, nrow = n.ahead, ncol = N^2)
  for (i in 2:(n.ahead+1)){
    tm2 <- sqrt(solve(diag(diag(H_t_upper[[i]]))))%*%H_t_upper[[i]]%*%sqrt(solve(diag(diag(H_t_upper[[i]]))))
    diag(tm2) <- sqrt(diag(H_t_upper[[i]]))
    sigma_t_upper[i-1,] <- c(tm2)
  }


  colnames(sigma_t_lower) <- rep(1, N^2)
  colnames(sigma_t_upper) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t_lower)[k2] <- paste('Lower CI conditional standard deviation of \n', colnames(x$data)[k])
        colnames(sigma_t_upper)[k2] <- paste('Upper CI conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t_lower)[k2] <- paste('Lower CI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        colnames(sigma_t_upper)[k2] <- paste('Upper CI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])

        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t_lower <- sigma_t_lower[, which(colSums(elim) == 1)]
  sigma_t_upper <- sigma_t_upper[, which(colSums(elim) == 1)]

  H_t_f_lower <- H_t_f_upper <- matrix(NA, nrow = n.ahead + 1, ncol = ncol(x$data)^2)

  for (i in 1:(n.ahead+1)){
    H_t_f_lower[i, ] <- c(H_t_lower[[i]])
    H_t_f_upper[i, ] <- c(H_t_upper[[i]])
  }

  if (inherits(x$data, "ts")) {
    sigma_t <- ts(sigma_t, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
    sigma_t_lower <- ts(sigma_t_lower, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
    sigma_t_upper <- ts(sigma_t_upper, start = (time(x$data)[nrow(x$data)]+1), frequency = frequency(x$data))
  }
  else if(inherits(x$data, "xts") || inherits(x$data, "zoo") ){
    sigma_t <- xts(matrix(sigma_t, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                   by = periodicity(x$data)$units))
    sigma_t_lower <- xts(matrix(sigma_t_lower, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                               by = periodicity(x$data)$units))
    sigma_t_upper <- xts(matrix(sigma_t_upper, nrow = n.ahead), order.by = seq((time(x$data)[nrow(x$data)]+1), (time(x$data)[nrow(x$data)] +  n.ahead),
                                                                               by = periodicity(x$data)$units))
  }


  result <- list(
    volatility_forecast = sigma_t,
    volatility_lower_conf_band = sigma_t_lower,
    volatility_upper_conf_band = sigma_t_upper,
    H_t_forecast = H_t_f,
    H_t_lower_conf_band = H_t_f_lower,
    H_t_upper_conf_band = H_t_f_upper,
    n.ahead = n.ahead,
    bekkfit = x
  )
  class(result) <-  c('bekkForecast', 'bekka')
  return(result)
}
