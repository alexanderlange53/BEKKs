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

#' @rdname bekk_fit_methods
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

#' @rdname bekk_fit_methods
#' @export
AIC.bekkFit <- function(x, ..., k = 2) {

  AICinner <- function(e) {
    N <- ncol(e$data)
    if (any(class(e) == 'bekk')) {
      aic <- k * 2 * N^2 + N * (N + 1)/2 - 2 * logLik(e)
    } else if (any(class(e) == 'bekka')) {
      aic <- k * 3 * N^2 + N * (N + 1)/2 - 2 * logLik(e)
    } else if (any(class(e) == 'sbekk')) {
      aic <- k * 2  + N * (N + 1)/2 - 2 * logLik(e)
    } else if (any(class(e) == 'sbekka')) {
      aic <- k * 3 + N * (N + 1)/2 - 2 * logLik(e)
    }else if (any(class(e) == 'dbekk')) {
      aic <- k * 2  + N * (N + 1)/2 - 2 * logLik(e)
    } else if (any(class(e) == 'dbekka')) {
      aic <- k * 3 + N * (N + 1)/2 - 2 * logLik(e)
    }
    return(aic)
  }

  if(!missing(...)) {# several objects: produce data.frame
    lls <- sapply(list(x, ...), AICinner)
    vals <- sapply(list(x, ...), function(e1){length(e1$theta)})
    aic <- data.frame(df = vals, AIC = lls)
  } else {
    aic <- AICinner(x)
  }

  return(aic)
}

#' @rdname bekk_fit_methods
#' @export
BIC.bekkFit <- function(x, ...) {

  BICinner <- function(e) {
    N <- ncol(e$data)
    if (any(class(e) == 'bekk')) {
      bic <- N^2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    } else if (any(class(e) == 'bekka')) {
      bic <- 3 * N^2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    } else if (any(class(e) == 'sbekk')) {
      bic <- 2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    } else if (any(class(e) == 'sbekka')) {
      bic <- 3 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    } else if (any(class(e) == 'dbekk')) {
      bic <- 2 * N + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    } else if (any(class(e) == 'dbekka')) {
      bic <- 3 * N + N * (N + 1)/2 * log(nrow(e$data)) - 2 * logLik(e)
    }
    return(bic)
  }

  if(!missing(...)) {# several objects: produce data.frame
    lls <- sapply(list(x, ...), BICinner)
    vals <- sapply(list(x, ...), function(e1){length(e1$theta)})
    bic <- data.frame(df = vals, BIC = lls)
  } else {
    bic <- BICinner(x)
  }

  return(bic)
}

#' @rdname bekk_fit_methods
#' @export
print.bekkFit <- function(x,...){
  bekkObject <- x

  if (any(class(bekkObject) == 'bekk')) {
    cat(paste("\n", "BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("BEKK estimation results")), collapse = "")
  } else if (any(class(bekkObject) == 'bekka')) {
    cat(paste("\n", "Asymmetric BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric BEKK estimation results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekk')) {
    cat(paste("\n", "Diagonal BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Diagonal BEKK estimation results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekka')) {
    cat(paste("\n", "Asymmetric diagonal BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric diagonal BEKK estimation results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekk')) {
    cat(paste("\n", "Scalar BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Scalar BEKK estimation results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekka')) {
    cat(paste("\n", "Asymmetric scalar BEKK estimation results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric scalar BEKK estimation results")), collapse = "")
  }

  cat(underScore)
  cat("\nLog-likelihood: ")
  cat(bekkObject$log_likelihood)
  cat("\nBEKK model stationary: ")
  cat(bekkObject$BEKK_valid)
  cat("\nNumber of BHHH iterations: ")
  cat(bekkObject$iter)
  cat("\nAIC: ")
  cat(bekkObject$AIC)
  cat("\nBIC: ")
  cat(bekkObject$BIC)
}

#' @rdname bekk_fit_methods
#' @export
residuals.bekkFit <- function(x, ...) {
  x$e_t
}
