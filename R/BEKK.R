#' Estimating a BEKK(1, 1) model
#'
#' @param r data input
#' @param init_values initial values for BEKK parameter
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#' @export

bekk <- function(r, init_values = NULL, max_iter = 300000) {
  N <- ncol(r)

  if (is.null(init_values)) {
      # Grid search
  } else {
      theta <- init_values
  }




}
