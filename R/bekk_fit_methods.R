#' bekkFit method
#'
#' @description Generic 'bekkFit' methods. More details on 'bekkFit' are described in \link{bekk_fit}
#'
#' @param x An object of class "bekkFit" from function \link{bekk_fit}.
#' @param object An object of class "bekkFit" from function \link{bekk_fit}.
#' @param k Numeric value, the penalty per parameter for AIC to be used; the default k = 2 is the classical AIC.
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
#' logLik(x1)
#' }
#' @import xts
#' @import stats
#' @rdname bekk_fit_methods
#' @export
print.bekkFit <- function(x,...){
  bekkObject <- x

  if (all(class(bekkObject) != 'var')) {
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
    cat("\n")
  }
  bekkObject
}

#' @rdname bekk_fit_methods
#' @export
residuals.bekkFit <- function(object, ...) {
  object$e_t
}




#' @noRd
AIC.bekkFit <- function(object, ..., k = 2) {

  x <- object

  AICinner <- function(e) {
    N <- ncol(e$data)
    if (any(class(e) == 'bekk')) {
      aic <- k * 2 * N^2 + N * (N + 1)/2 - 2 * llv(e)
    } else if (any(class(e) == 'bekka')) {
      aic <- k * 3 * N^2 + N * (N + 1)/2 - 2 * llv(e)
    } else if (any(class(e) == 'sbekk')) {
      aic <- k * 2  + N * (N + 1)/2 - 2 * llv(e)
    } else if (any(class(e) == 'sbekka')) {
      aic <- k * 3 + N * (N + 1)/2 - 2 * llv(e)
    }else if (any(class(e) == 'dbekk')) {
      aic <- k * 2  + N * (N + 1)/2 - 2 * llv(e)
    } else if (any(class(e) == 'dbekka')) {
      aic <- k * 3 + N * (N + 1)/2 - 2 * llv(e)
    }
    return(aic)
  }

  if(!missing(...)) {# several objects: produce data.frame
    lls <- sapply(list(x, ...), AICinner)
    vals <- sapply(list(x, ...), function(e1){length(e1$theta)})
    aic <- data.frame(df = vals, AIC = lls)
  } else {
    aic <- AICinner(x)
    aic <- data.frame(df = length(x$theta), AIC = aic)
  }

  return(aic)
}
#' @noRd
BIC.bekkFit <- function(object, ...) {

  x <- object

  BICinner <- function(e) {
    N <- ncol(e$data)
    if (any(class(e) == 'bekk')) {
      bic <- 2* N^2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    } else if (any(class(e) == 'bekka')) {
      bic <- 3 * N^2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    } else if (any(class(e) == 'sbekk')) {
      bic <- 2 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    } else if (any(class(e) == 'sbekka')) {
      bic <- 3 + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    } else if (any(class(e) == 'dbekk')) {
      bic <- 2 * N + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    } else if (any(class(e) == 'dbekka')) {
      bic <- 3 * N + N * (N + 1)/2 * log(nrow(e$data)) - 2 * llv(e)
    }
    return(bic)
  }

  if(!missing(...)) {# several objects: produce data.frame
    lls <- sapply(list(x, ...), BICinner)
    vals <- sapply(list(x, ...), function(e1){length(e1$theta)})
    bic <- data.frame(df = vals, BIC = lls)
  } else {
    bic <- data.frame(df = length(x$theta), BIC = BICinner(x))
  }

  return(bic)
}

#' @noRd
llv <- function(object) {


  llv_inner <- function(x){
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


    lls <- llv_inner(object)



  return(lls)
}

#' @rdname bekk_fit_methods
#' @export
logLik.bekkFit <- function(object, ..., k = 2) {


  llv_inner <- function(x){
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

  if(!missing(...)) {# several objects: produce data.frame
    lls <- sapply(list(object, ...), llv_inner)
    vals <- sapply(list(object, ...), function(e1){length(e1$theta)})
    aic <- data.frame(df = vals, LLV = lls, AIC = AIC(object, ..., k = k)$AIC, BIC = BIC(object, ...)$BIC)
  } else {
    lls <- llv_inner(object)
    aic <- data.frame(df = length(object$theta), LLV = lls, AIC = AIC(object, k = k)$AIC, BIC = BIC(object)$BIC)


     }
  return(aic)
}

