#' Simulating BEKK models
#'
#' @description Method for simulating a N-dimensional BEKK model.
#'
#' @param spec A spec object of class "bekkSpec" from the function \link{bekk_spec} or a fitted bekk model of class "bekkFit" from the \link{bekk_fit} function
#' @param nobs Number of observations of the simulated sample
#'
#'
#' @examples
#' \donttest{
#'
#' # Simulate a BEKK with estimated parameter
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds)
#'
#' x2 <- bekk_sim(x1, 3000)
#'
#' plot(x2)
#'
#' }
#'
#' @export

bekk_sim <- function(spec, nobs) {
  UseMethod('bekk_sim')
}

#' @export
bekk_sim.bekk <- function(spec, nobs) {

  xx <- process_object(spec)
  par <- coef_mat(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    cat("Please provide a stationary BEKK model")
    return(NULL)
  }

  sim_dat <- simulate_bekk_c(c(xx$theta), nobs, xx$N)

  return(ts(sim_dat))
}

#' @export
bekk_sim.bekka <- function(spec, nobs) {

  xx <- process_object(spec)
  par <- coef_mat_asymm(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    cat("Please provide a stationary BEKK model")
    return(NULL)
  }
  #expected_signs
  sim_dat <- simulate_bekka_c(c(xx$theta), nobs, xx$N, xx$signs, xx$expected_signs)

  return(ts(sim_dat))
}

