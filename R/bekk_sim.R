#' Simulating BEKK models
#'
#' @param spec A spec object of class bekkSimSpec or a fitted bekk model of class bekk from the \link{bekk} function
#' @param nobs Number of observations of the simulated sample
#'
#'
#' @export

bekk_sim <- function(spec, nobs) {
  UseMethod('bekk_sim')
}

#' @export
bekk_sim.bekk <- function(spec, nobs) {

  theta <- spec$theta

  N <- ncol(spec$C0)

  sim_dat <- simulate_bekk_c(c(theta), nobs, N)


  return(ts(sim_dat))
}

# if (inherits(spec, 'bekk')) {
#   theta <- spec$theta
#   N <- ncol(spec$C0)
#
#   sim_dat <- simulate_bekk_c(c(theta), nobs, N)
# } else if (class(spec) == 'bekkSimSpec') {
#
# } else {
#   stop("The object 'spec' should be of class 'bekk' or 'bekkSimSpec'.")
# }
