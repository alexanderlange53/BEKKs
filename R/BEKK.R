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
#' }
#' @export

bekk <- function(r, init_values = NULL,
                 QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9){
  N <- ncol(r)

  r <- data.matrix(r)

  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_BEKK(r)
    } else if (init_values == 'random') {

    }
  } else {
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
  elim <- elimination_mat(N)
  sigma_t <- t(elim %*% t(var_process$sigma_t))

  return(list(C0 =  param_mat$c0,
              A = param_mat$a,
              G = param_mat$g,
              C0_t = tratios_mat$c0,
              A_t = tratios_mat$a,
              G_t = tratios_mat$g,
              log_likelihood = params$likelihood,
              sigma_t = sigma_t,
              e_t = var_process$e_t,
              iter = params$iter))
}