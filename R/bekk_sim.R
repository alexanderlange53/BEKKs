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

  xx <- process_object(spec)

  sim_dat <- simulate_bekk_c(c(xx$theta), nobs, xx$N)

  return(ts(sim_dat))
}

#' @export
bekk_sim.bekka <- function(spec, nobs) {

  xx <- process_object(spec)

  sim_dat <- simulate_bekka_c(c(xx$theta), nobs, xx$N, xx$signs)

  return(ts(sim_dat))
}

