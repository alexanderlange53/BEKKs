#' @name simulate
#' @rdname simulate
#' @title Simulating BEKK models
#'
#' @description Method for simulating a N-dimensional BEKK model.
#'
#' @param object A spec object of class "bekkSpec" from the function \link{bekk_spec} or a fitted bekk model of class "bekkFit" from the \link{bekk_fit} function
#' @param nsim Number of observations of the simulated sample
#' @param ... Further parameters to be passed on to the function.
#' @return Returns a simulated time series S3 class object using the parameters of passed "bekkSpec" or "bekkFit".
#'
#' @examples
#' \donttest{
#'
#' # simulate a BEKK with estimated parameter
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds)
#'
#' x2 <- simulate(x1, nsim = 3000)
#'
#' plot(x2)
#'
#' }
#'

#' @export
simulate.bekk <- function(object, nsim, ...) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")

  }

  sim_dat <- simulate_bekk_c(c(xx$theta), nsim, xx$N)

  return(ts(sim_dat))
}
#' @rdname simulate

#' @export
simulate.bekka <- function(object, ..., nsim) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat_asymm(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")

  }
  #expected_signs
  sim_dat <- simulate_bekka_c(c(xx$theta), nsim, xx$N, xx$signs, xx$expected_signs)

  return(ts(sim_dat))
}

#' @rdname simulate

#' @export
simulate.dbekk <- function(object, ..., nsim) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat_diagonal(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")

  }

  sim_dat <- simulate_dbekk_c(c(xx$theta), nsim, xx$N)

  return(ts(sim_dat))
}
#' @rdname simulate
#' @export
simulate.dbekka <- function(object, ..., nsim) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat_asymm_diagonal(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")
  }
  #expected_signs
  sim_dat <- simulate_dbekka_c(c(xx$theta), nsim, xx$N, xx$signs, xx$expected_signs)

  return(ts(sim_dat))
}

#' @rdname simulate
#' @export
simulate.sbekk <- function(object, ..., nsim) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat_scalar(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")
  }

  sim_dat <- simulate_sbekk_c(c(xx$theta), nsim, xx$N)

  return(ts(sim_dat))
}
#' @rdname simulate
#' @export
simulate.sbekka <- function(object, ..., nsim) {
  spec <- object
  if(is.null(nsim) || !is.numeric(nsim) || nsim < 1){
    stop("Please provide an integer specifying the number of observations")
  }
  if(!inherits(spec,c("bekkSpec", "bekkFit"))){
    stop("Please provide an object of class bekk_fit or 'bekk_spec'.")
  }
  xx <- process_object(spec)
  par <- coef_mat_asymm_scalar(xx$theta,xx$N)
  if(xx$BEKK_valid==FALSE){
    stop("Please provide a stationary BEKK model.")
      }
  #expected_signs
  sim_dat <- simulate_sbekka_c(c(xx$theta), nsim, xx$N, xx$signs, xx$expected_signs)

  return(ts(sim_dat))
}

