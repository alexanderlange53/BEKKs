#' Estimating a BEKK(1, 1) model
#'
#' @param r data input
#' @param init_values initial values for BEKK parameter
#' @param QML_t_ratios Logical. If QML_t_ratios = 'TRUE', the t-ratios of the BEKK parameter matrices
#'                     are exactly calculated via second order derivatives.
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#'
#' @examples
#' \donttest{
#'
#' data(bivariate)
#' x1 <- bekk(BI, init_values = NULL,
#' QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' }
#' @export

bekk <- function(r, init_values = NULL,
                 QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9){

  # Checking for valid input
  r <- as.matrix(r)
  if (any(is.na(r)))
    stop("\nNAs in y.\n")
  if (ncol(r) < 2)
    stop("The matrix 'r' should contain at least two variables.")
  if (is.null(colnames(r))) {
    colnames(r) <- paste("y", 1:ncol(r), sep = "")
  }


  N <- ncol(r)

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_BEKK(r)
    } else if (init_values == 'random') {
      theta <- random_grid_search_BEKK(r, 1000)
      theta <- theta[[1]]
    }
  } else {
    if(length(init_values) != 2 * N^2 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_bekk(r, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios(params$theta, r)
    tratios_mat <- coef_mat(tratios, N)
  } else {
    tratios_mat <- coef_mat(params$t_val, N)
  }

  param_mat <- coef_mat(params$theta, N)

  var_process <- sigma_bekk(r, param_mat$c0, param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional variances of', colnames(r)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional covariances of', colnames(r)[i], ' and ', colnames(r)[j])
        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  # Final check if BEKK is valid
  BEKK_valid <- valid_bekk(param_mat$c0, param_mat$a, param_mat$g)


  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 G_t = tratios_mat$g,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 e_t = var_process$e_t,
                 iter = params$iter,
                 data = r)
  class(result) <- 'bekk'
  return(result)
}
