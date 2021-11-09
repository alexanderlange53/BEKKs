#' Estimating BEKK(1, 1) models
#'
#' @param spec object of class "bekkSpec" from function \href{bekk_spec}.
#' @param data data input
#' @param QML_t_ratios Logical. If QML_t_ratios = 'TRUE', the t-ratios of the BEKK parameter matrices
#'                     are exactly calculated via second order derivatives.
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#'
#' @examples
#' \donttest{
#'
#' data(bivariate)
#' x1 <- bekk_fit(BI, init_values = NULL,
#' QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' }
#' @export

bekk_fit <- function(spec, data, QML_t_ratios = FALSE,
                 seed = NULL, max_iter = 50, crit = 1e-9, nc = 1){



  # Checking for valid input
  if (class(spec) != 'bekkSpec') {
    stop("Please provide an object of class bekkSpec.")
  }

  r <- as.matrix(data)
  if (any(is.na(r))) {
    stop("\nNAs in data.\n")
  }
  if (ncol(r) < 2) {
    stop("The data matrix should contain at least two variables.")
  }
  if (is.null(colnames(r))) {
    colnames(r) <- paste("y", 1:ncol(r), sep = "")
  }


  N <- ncol(r)

  if(spec[[1]]$type == "bekk") {
    if(spec[[1]]$asymmetric == FALSE){
      UseMethod("bekk")
    }else {
      UseMethod("bekka")
    }
  }


}
