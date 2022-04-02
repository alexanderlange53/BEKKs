#' Estimating multivariate BEKK-type volatility models
#'
#' @description Method for fitting a variety of N-dimensional BEKK models.
#'
#' @param spec An object of class "bekkSpec" from function \link{bekk_spec}.
#' @param data A multivariate data object. Can be a numeric matrix or ts/xts/zoo object.
#' @param QML_t_ratios Logical. If QML_t_ratios = 'TRUE', the t-ratios of the BEKK parameter matrices
#'                     are exactly calculated via second order derivatives.
#' @param max_iter Maximum number of BHHH algorithm iterations.
#' @param crit Determines the precision of the BHHH algorithm.
#' @return  Returns a S3 class "bekkFit" object containing the estimated parameters, t-values, volatility process of the model defined by the BEKK_spec object.
#'
#' @details The BEKK optimization routine is based on the Berndt–Hall–Hall–Hausman (BHHH) algorithm and is inspired by the study of Hafner and Herwartz (2008).
#' The authors provide analytical formulas for the score and Hessian of several MGARCH models in a QML framework and show that analytical derivations significantly outperform numerical methods.
#'
#' @references Hafner and Herwartz (2008). Analytical quasi maximum likelihood inference in multivariate volatility models. Metrika, 67, 219-239.
#'
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#'
#' # Fitting a symmetric BEKK model
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' # Fitting an asymmetric BEKK model
#' obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#' x1 <- bekk_fit(obj_spec, StocksBonds)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' # Fitting a symmetric diagonal BEKK model
#' obj_spec <- bekk_spec(model = list(type = "dbekk", asymmetric = FALSE))
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#'
#' # Fitting a symmetric scalar BEKK model
#' obj_spec <- bekk_spec(model = list(type = "sbekk", asymmetric = FALSE))
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' }
#' @import xts
#' @import stats
#' @export

bekk_fit <- function(spec, data, QML_t_ratios = FALSE,
                     max_iter = 50, crit = 1e-9){

  if (!inherits(spec, 'bekkSpec')) {
    stop('Please provide and object of class "bekkSpec" for spec.')
  }

  if (any(is.na(data))) {
    stop("\nNAs in data.\n")
  }
  if (ncol(data) < 2) {
    stop("The data matrix should contain at least two variables.")
  }
  if (is.null(colnames(data))) {
    colnames(data) <- paste("y", 1:ncol(data), sep = "")
  }

  UseMethod('bekk_fit')

}

#' @export
bekk_fit.bekk <- function(spec, data, QML_t_ratios = FALSE,
                          max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_BEKK(data)
      theta <- theta[[1]]
    } else if (init_values == 'random') {
      cat('Generating starting values \n')

      theta <- random_grid_search_BEKK(data)
      theta <- theta[[1]]
    } else if (init_values == 'simple') {
      uncond_var <- crossprod(data)/nrow(data)
      A <- matrix(0, ncol = N, nrow = N)
      G <- matrix(0, ncol = N, nrow = N)
      C <- matrix(0, ncol = N, nrow = N)
      #th0=numeric(2*n^2+n*(n+1)/2)

      diag(A) <- 0.3
      diag(G) <- 0.92
      diag(C) <- 0.05*diag(uncond_var)


      for (i in 1:N){
        for (j in seq(i,N)){

          cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
          C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
          C[j,i] <- C[i, j]

        }
      }

      C = t(chol(C))
      C0 = C[,1]

      if (N > 2) {
        for (i in 2:(N-1)){
          C0 = c(C0, C[i:N, i])
        }
      }

      C0 = c(C0, C[N, N])

      theta = c(C0, c(A), c(G))
     }
  } else {
    if(length(init_values) != 2 * N^2 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_bekk(data, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios(params$theta, data)
    tratios_mat <- coef_mat(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat(abs(params$t_val), N)
  }

  param_mat <- coef_mat(params$theta, N)


  var_process <- sigma_bekk(data, t(param_mat$c0), param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]


  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_bekk(param_mat$c0, param_mat$a, param_mat$g)


  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 G_t = tratios_mat$g,
                 theta = params$theta,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = FALSE,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'bekk')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}

#' @export
bekk_fit.bekka <- function(spec, data, QML_t_ratios = FALSE,
                   max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(is.null(spec$model$signs)){
    spec$model$signs = matrix(rep(-1, N), ncol = 1)
  }
  if(length(spec$model$signs) != N){
    stop('Length of "signs" does not match dimension of data.')
  }

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_asymmetricBEKK(data, spec$model$signs)
      theta <- theta[[1]]
    } else if (init_values == 'random') {

      cat('Generating starting values \n')
      theta = random_grid_search_asymmetric_BEKK(data, spec$model$signs)[[1]]
    } else if (init_values == 'simple') {
      uncond_var <- crossprod(data)/nrow(data)
      A <- matrix(0, ncol = N, nrow = N)
      B <- matrix(0, ncol = N, nrow = N)
      G <- matrix(0, ncol = N, nrow = N)
      C <- matrix(0, ncol = N, nrow = N)
      #th0=numeric(2*n^2+n*(n+1)/2)

      diag(A) <- 0.25
      diag(B) <- 0.05
      diag(G) <- 0.92
      diag(C) <- 0.05*diag(uncond_var)


      for (i in 1:N){
        for (j in seq(i,N)){

          cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
          C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
          C[j,i] <- C[i, j]

        }
      }

      C = t(chol(C))
      C0 = C[,1]

      if (N > 2) {
        for (i in 2:(N-1)){
          C0 = c(C0, C[i:N, i])
        }
      }

      C0 = c(C0, C[N, N])

      theta = c(C0, c(A), c(G))

    }
  } else {
    if(length(init_values) != 3 * N^2 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_asymm_bekk(data, theta, max_iter, crit, spec$model$signs)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios_asymm(params$theta, data, spec$model$signs)
    tratios_mat <- coef_mat_asymm(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat_asymm(abs(params$t_val), N)
  }

  param_mat <- coef_mat_asymm(params$theta, N)

  var_process <- sigma_bekk_asymm(data, t(param_mat$c0), param_mat$a, param_mat$b, param_mat$g, spec$model$signs)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional SD of', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_asymm_bekk(param_mat$c0, param_mat$a, param_mat$b, param_mat$g, data, spec$model$signs)

  expected_signs=expected_indicator_value(data,spec$model$signs)

  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 B = param_mat$b,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 B_t = tratios_mat$b,
                 G_t = tratios_mat$g,
                 theta = params$theta,
                 signs = spec$model$signs,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = TRUE,
                 expected_signs = expected_signs,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'bekka')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}

#' @export
bekk_fit.dbekk <- function(spec, data, QML_t_ratios = FALSE,
                          max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_dBEKK(data)
      theta <- theta[[1]]
    } else if (init_values == 'random') {
      cat('Generating starting values \n')

      theta <- random_grid_search_dBEKK(data)
      theta <- theta[[1]]
    }
  } else {
    if(length(init_values) != 2 * N + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_dbekk(data, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios_dbekk(params$theta, data)
    tratios_mat <- coef_mat_diagonal(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat_diagonal(abs(params$t_val), N)
  }

  param_mat <- coef_mat_diagonal(params$theta, N)


  var_process <- sigma_bekk(data, t(param_mat$c0), param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]


  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_bekk(param_mat$c0, param_mat$a, param_mat$g)


  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 G_t = tratios_mat$g,
                 theta = params$theta,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = FALSE,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'dbekk')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}

#' @export
bekk_fit.dbekka <- function(spec, data, QML_t_ratios = FALSE,
                           max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(is.null(spec$model$signs)){
    spec$model$signs = matrix(rep(-1, N), ncol = 1)
  }
  if(length(spec$model$signs) != N){
    stop('Length of "signs" does not match dimension of data.')
  }

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_asymmetricdBEKK(data, spec$model$signs)
      theta <- theta[[1]]
    } else if (init_values == 'random') {

      cat('Generating starting values \n')
      theta = random_grid_search_asymmetric_dBEKK(data, spec$model$signs)[[1]]
    }
  } else {
    if(length(init_values) != 3 * N + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_asymm_dbekk(data, theta, max_iter, crit, spec$model$signs)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios_dbekka(params$theta, data, spec$model$signs)
    tratios_mat <- coef_mat_asymm_diagonal(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat_asymm_diagonal(abs(params$t_val), N)
  }

  param_mat <- coef_mat_asymm_diagonal(params$theta, N)

  var_process <- sigma_bekk_asymm(data, t(param_mat$c0), param_mat$a, param_mat$b, param_mat$g, spec$model$signs)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional SD of', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_asymm_bekk(param_mat$c0, param_mat$a, param_mat$b, param_mat$g, data, spec$model$signs)

  expected_signs=expected_indicator_value(data,spec$model$signs)

  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 B = param_mat$b,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 B_t = tratios_mat$b,
                 G_t = tratios_mat$g,
                 theta = params$theta,
                 signs = spec$model$signs,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = TRUE,
                 expected_signs = expected_signs,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'dbekka')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}

#' @export
bekk_fit.sbekk <- function(spec, data, QML_t_ratios = FALSE,
                          max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_sBEKK(data)
      theta <- theta[[1]]
    } else if (init_values == 'random') {
      cat('Generating starting values \n')

      theta <- random_grid_search_sBEKK(data)
      theta <- theta[[1]]
    }
  } else {
    if(length(init_values) != 2  + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_sbekk(data, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios_sbekk(params$theta, data)
    tratios_mat <- coef_mat_scalar(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat_scalar(abs(params$t_val), N)
  }

  param_mat <- coef_mat_scalar(params$theta, N)


  var_process <- sigma_sbekk(data, t(param_mat$c0), param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]


  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_sbekk(param_mat$c0, param_mat$a, param_mat$g)


  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 a = param_mat$a,
                 g = param_mat$g,
                 C0_t = tratios_mat$c0,
                 a_t = tratios_mat$a,
                 g_t = tratios_mat$g,
                 theta = params$theta,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = FALSE,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'sbekk')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}

#' @export
bekk_fit.sbekka <- function(spec, data, QML_t_ratios = FALSE,
                           max_iter = 50, crit = 1e-9) {

  init_values <- spec$init_values
  N <- ncol(data)

  if(is.null(spec$model$signs)){
    spec$model$signs = matrix(rep(-1, N), ncol = 1)
  }
  if(length(spec$model$signs) != N){
    stop('Length of "signs" does not match dimension of data.')
  }

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_asymmetricsBEKK(data, spec$model$signs)
      theta <- theta[[1]]
    } else if (init_values == 'random') {

      cat('Generating starting values \n')
      theta = random_grid_search_asymmetric_sBEKK(data, spec$model$signs)[[1]]
    }
  } else {
    if(length(init_values) != 3 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_asymm_sbekk(data, theta, max_iter, crit, spec$model$signs)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios_sbekk_asymm(params$theta, data, spec$model$signs)
    tratios_mat <- coef_mat_asymm_scalar(abs(tratios), N)
  } else {
    tratios_mat <- coef_mat_asymm_scalar(abs(params$t_val), N)
  }

  param_mat <- coef_mat_asymm_scalar(params$theta, N)

  var_process <- sigma_sbekk_asymm(data, t(param_mat$c0), param_mat$a, param_mat$b, param_mat$g, spec$model$signs)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional SD of', colnames(data)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of', colnames(data)[i], ' and ', colnames(data)[j])
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  if (inherits(data, "ts")) {
    sigma_t <- ts(sigma_t, start = time(data)[1], frequency = frequency(data))
  }
  else if(inherits(data, "xts") || inherits(data, "zoo") ){
    sigma_t <- xts(sigma_t, order.by = time(data))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_asymm_sbekk(param_mat$c0, param_mat$a, param_mat$b, param_mat$g, data, spec$model$signs)

  expected_signs=expected_indicator_value(data,spec$model$signs)

  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 a = param_mat$a,
                 b = param_mat$b,
                 g = param_mat$g,
                 C0_t = tratios_mat$c0,
                 a_t = tratios_mat$a,
                 b_t = tratios_mat$b,
                 g_t = tratios_mat$g,
                 theta = params$theta,
                 signs = spec$model$signs,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 H_t = var_process$sigma_t,
                 e_t = var_process$e_t,
                 Second_moments_of_residuals = cov(var_process$e_t),
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 asymmetric = TRUE,
                 expected_signs = expected_signs,
                 data = data,
                 spec = spec)
  class(result) <- c('bekkFit', 'sbekka')

  result$AIC <- AIC(result)
  result$BIC <- BIC(result)

  return(result)
}
