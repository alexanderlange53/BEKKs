#' Forecasting conditional volatilities with BEKK models
#'
#' @param x A fitted bekk model of class bekk from the \link{bekk} function
#' @param n.ahead Number of periods to forecast conditional volatility. Default is a one-period ahead forecast.
#' @param KI_niveau Floating point in [0,1] defining the niveau for confidence bands of the conditional volatility forecast. Provided are either 90%, 95% or 99% confidence levels. Default are 95% niveau confidence bands.
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
bekk_forecast <- function(x, n.ahead = 1, KI_niveau=0.95){

  if (!inherits(x, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for "x".')
  }


  UseMethod('bekk_forecast')

}

#' @export
bekk_forecast.bekk <- function(x, n.ahead = 1, KI_niveau=0.95) {
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

  #95% confidence interval

  score_final = score_bekk(x$theta, x$data)
  s1_temp = diag(solve(t(score_final) %*% score_final),names=T)
  s1 = sqrt(s1_temp)

  if(KI_niveau == 0.95){
  lower_theta = x$theta-1.64*s1
  upper_theta = x$theta+1.64*s1
  }
  else if(KI_niveau == 0.99){
    lower_theta = x$theta-2.58*s1
    upper_theta = x$theta+2.58*s1
  } else{
    lower_theta = x$theta-1.96*s1
    upper_theta = x$theta+1.96*s1
  }

  H_t_lower <- vector(mode = "list",length = n.ahead+1)
  H_t_upper <- vector(mode = "list",length = n.ahead+1)
  x_lower <- coef_mat(lower_theta, N)
  x_upper <- coef_mat(upper_theta, N)
  #check if t(x_lower$C0) or x_lower$C0 for sigma_bekk
  H_t_lower[[1]] <- matrix(sigma_bekk(x$data,t(x_lower$C0),x_lower$A,x_lower$G)$sigma_t[NoBs,],nrow = N, ncol = N)
  H_t_upper[[1]] <- matrix(sigma_bekk(x$data,t(x_upper$C0),x_upper$A,x_upper$G)$sigma_t[NoBs,],nrow = N, ncol = N)
  current_returns <- t(x$data[NoBs,])


  for(i in 1:n.ahead){
    H_t_lower[[i+1]] <- t(x_lower$C0) %*% x_lower$C0 + t(x_lower$A) %*% t(current_returns) %*% current_returns %*% x_lower$A + t(x_lower$G) %*% H_t[[i]] %*% x_lower$G
    current_returns <- eigen_value_decomposition(H_t_lower[[i+1]])
  }

  current_returns <- t(x$data[NoBs,])

  for(i in 1:n.ahead){
    H_t_upper[[i+1]] <- t(x_upper$C0) %*% x_upper$C0 + t(x_upper$A) %*% t(current_returns) %*% current_returns %*% x_upper$A + t(x_upper$G) %*% H_t[[i]] %*% x_upper$G
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
        colnames(sigma_t_lower)[k2] <- paste('Lower KI conditional standard deviation of \n', colnames(x$data)[k])
        colnames(sigma_t_upper)[k2] <- paste('Upper KI conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t_lower)[k2] <- paste('Lower KI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        colnames(sigma_t_upper)[k2] <- paste('Upper KI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])

        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t_lower <- sigma_t_lower[, which(colSums(elim) == 1)]
  sigma_t_upper <- sigma_t_upper[, which(colSums(elim) == 1)]


  result <- list(
    volatility_forecast = sigma_t,
    H_t_forecast = H_t_f,
    volatility_lower_conf_band = sigma_t_lower,
    volatility_upper_conf_band = sigma_t_upper,
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

  expected_signs=expected_indicator_value(x$data,x$signs)
  for(i in 2:n.ahead){
    H_t[[i+1]]=t(x$C0) %*% x$C0 + t(x$A) %*% H_t[[i+1]] %*% x$A + expected_signs * t(x$B) %*% H_t[[i]] %*% x$B + t(x$G) %*% H_t[[i]] %*% x$G

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

  if(KI_niveau == 0.95){
    lower_theta = x$theta-1.64*s1
    upper_theta = x$theta+1.64*s1
  }
  else if(KI_niveau == 0.99){
    lower_theta = x$theta-2.58*s1
    upper_theta = x$theta+2.58*s1
  } else{
    lower_theta = x$theta-1.96*s1
    upper_theta = x$theta+1.96*s1
  }

  H_t_lower <- vector(mode = "list",length = n.ahead+1)
  H_t_upper <- vector(mode = "list",length = n.ahead+1)
  x_lower <- coef_mat_asymm(lower_theta, N)
  x_upper <- coef_mat_asymm(upper_theta, N)
  #check if t(x_lower$C0) or x_lower$C0 for sigma_bekk
  H_t_lower[[1]] <- matrix(sigma_bekk_asymm(x$data, t(x_lower$C0), x_lower$A, x_lower$B, x_lower$G, x$signs)$sigma_t[NoBs,],nrow = N, ncol = N)
  H_t_upper[[1]] <- matrix(sigma_bekk_asymm(x$data, t(x_upper$C0), x_upper$A, x_upper$B, x_upper$G, x$signs)$sigma_t[NoBs,],nrow = N, ncol = N)
  current_returns <- t(x$data[NoBs,])


  for(i in 1:n.ahead){
    H_t_lower[[i+1]] <- t(x_lower$C0) %*% x_lower$C0 + t(x_lower$A) %*% t(current_returns) %*% current_returns %*% x_lower$A + expected_signs * t(x_lower$B) %*% H_t_lower[[i]] %*% x_lower$B + t(x_lower$G) %*% H_t[[i]] %*% x_lower$G
    current_returns <- eigen_value_decomposition(H_t_lower[[i+1]])
  }

  current_returns <- t(x$data[NoBs,])

  for(i in 1:n.ahead){
    H_t_upper[[i+1]] <- t(x_upper$C0) %*% x_upper$C0 + t(x_upper$A) %*% t(current_returns) %*% current_returns %*% x_upper$A + expected_signs * t(x_upper$B) %*% H_t_upper[[i]] %*% x_upper$B + t(x_upper$G) %*% H_t[[i]] %*% x_upper$G
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
        colnames(sigma_t_lower)[k2] <- paste('Lower KI conditional standard deviation of \n', colnames(x$data)[k])
        colnames(sigma_t_upper)[k2] <- paste('Upper KI conditional standard deviation of \n', colnames(x$data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t_lower)[k2] <- paste('Lower KI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])
        colnames(sigma_t_upper)[k2] <- paste('Upper KI conditional correlation of \n', colnames(x$data)[i], ' and ', colnames(x$data)[j])

        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t_lower <- sigma_t_lower[, which(colSums(elim) == 1)]
  sigma_t_upper <- sigma_t_upper[, which(colSums(elim) == 1)]


  result <- list(
    volatility_forecast = sigma_t,
    H_t_forecast = H_t_f,
    volatility_lower_conf_band = sigma_t_lower,
    volatility_upper_conf_band = sigma_t_upper,
    n.ahead = n.ahead,
    bekkfit = x
  )
  class(result) <-  c('bekkForecast', 'bekka')
  return(result)
}
