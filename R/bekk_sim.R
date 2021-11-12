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

#' @export
bekk_sim.bekkSpec <- function(spec, nobs) {

  theta <- spec$init_values

  N <- spec$N

  sim_dat <- simulate_bekk_c(c(theta), nobs, N)


  return(ts(sim_dat))
}
