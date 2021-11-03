#' Simulating BEKK models
#'
#' @param spec A spec object of class bekkSimSpec or a fitted bekk model of class bekk from the \link{bekk} function
#' @param nobs Number of observations of the simulated sample
#'
#'
#' @export

bekk_sim <- function(spec, nobs) {

  if (class(spec) == 'bekk') {
    theta <- spec$theta
    N <- ncol(spec$C0)

    sim_dat <- simulate_bekk_c(c(theta), nobs, N)
  } else if (class(spec) == 'bekkSimSpec') {

  } else {
    stop("The object 'spec' should be of class 'bekk' or 'bekkSimSpec'.")
  }

  return(ts(sim_dat))
}
