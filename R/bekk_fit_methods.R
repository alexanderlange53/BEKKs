#' bekkFit method
#'
#' @description Generic 'bekkFit' methods. More details on 'bekkFit' are described in \link{bekk_fit}
#'
#' @param x An object of class "bekkFit" from function \link{bekk_fit}.
#' @param k Numeric value, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ...	Further arguments to be passed to and from other methods.
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
#' AIC(x1)
#' }
#' @import xts
#' @import stats

#' @export
logLik.bekkFit <- function(x, ....) {

  if (any(class(x) == 'bekk')) {
    logl <- loglike_bekk(x$theta, x$data)
  } else if (any(class(x) == 'bekka')) {
    logl <- loglike_asymm_bekk(x$theta, x$data, x$signs)
  } else if (any(class(x) == 'sbekk')) {
    logl <- loglike_sbekk(x$theta, x$data)
  } else if (any(class(x) == 'sbekka')) {
    logl <- loglike_asymm_sbekk(x$theta, x$data, x$signs)
  }else if (any(class(x) == 'dbekk')) {
    logl <- loglike_dbekk(x$theta, x$data)
  } else if (any(class(x) == 'dbekka')) {
    logl <- loglike_asymm_dbekk(x$theta, x$data, x$signs)
  }

  return(logl)
}

#' @export
AIC.bekkFit <- function(x, ..., k = 2) {
  N <- ncol(x$data)
  if (any(class(x) == 'bekk')) {
    aic <- k * 2 * N^2 + N * (N + 1)/2 - 2 * logLik(x)
  } else if (any(class(x) == 'bekka')) {
    aic <- k * 3 * N^2 + N * (N + 1)/2 - 2 * logLik(x)
  } else if (any(class(x) == 'sbekk')) {
    aic <- k * 2  + N * (N + 1)/2 - 2 * logLik(x)
  } else if (any(class(x) == 'sbekka')) {
    aic <- k * 3 + N * (N + 1)/2 - 2 * logLik(x)
  }else if (any(class(x) == 'dbekk')) {
    aic <- k * 2  + N * (N + 1)/2 - 2 * logLik(x)
  } else if (any(class(x) == 'dbekka')) {
    aic <- k * 3 + N * (N + 1)/2 - 2 * logLik(x)
  }

  return(aic)
}

#' @export
BIC.bekkFit <- function(x, ...) {
  N <- ncol(x$data)
  if (any(class(x) == 'bekk')) {
    bic <- N^2 + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  } else if (any(class(x) == 'bekka')) {
    bic <- 3 * N^2 + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  } else if (any(class(x) == 'sbekk')) {
    bic <- 2 + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  } else if (any(class(x) == 'sbekka')) {
    bic <- 3 + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  } else if (any(class(x) == 'dbekk')) {
    bic <- 2 * N + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  } else if (any(class(x) == 'dbekka')) {
    bic <- 3 * N + N * (N + 1)/2 * log(nrow(x$data)) - 2 * logLik(x)
  }


  return(bic)
}
