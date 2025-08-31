#' BEKK specification method
#'
#' @description Method for creating a N-dimensional BEKK model specification object prior to fitting and/or simulating.
#'
#' @param model A list containing the model type specification: Either "bekk" "dbekk" or "sbekk".
#' Moreover it can be specified whether the model should be estimated allowing for asymmetric volatility structure.
#' @param init_values initial values for \link{bekk_fit} during BHHH algorithm. It can be either a numerical vector of suitable dimension, 'NULL' (default) to use a simple grid search algorithm, or a character vector i.e. "random" to use a random starting value generator (set a seed in advance for reproducible results), or
#'  "simple" for relying on a simple initial values generator based on typical values for BEKK parameter found in the literature. If the object from this function is passed to \link{simulate}, "init_values" are used as parameters for data generating process.
#' @param signs An N-dimensional vector consisting of "1" or "-1" to indicate the asymmetric effects to be considered.
#' Setting the i-th element of the vector to "1" or "-1" means that the model takes into account additional volatility if the returns of the i-th column in the data matrix are either positive or negative.
#' If "asymmetric = TRUE", the default is set to "rep(-1, N)" i.e. it is assumed that excess volatility occurs for all series if the returns are negative.
#' @param N Integer specifying the dimension of the BEKK model. Only relevant when this object of class "bekkSpec"" is used for simulating BEKK processes by applying it to \link{simulate}.
#' @return Returns a S3 class "bekkSpec"  object containing the specifications of the model to be estimated.
#'
#'
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#'
#' # Fitting a symmetric BEKK model using default starting values
#' # - i.e. fixed values
#' obj_spec_fixed <- bekk_spec(init_values = NULL)
#' x1 <- bekk_fit(obj_spec_fixed, StocksBonds, QML_t_ratios = FALSE,
#' max_iter = 50, crit = 1e-9)
#'# Fitting a symmetric BEKK model using initial values originating from a
#'# random grid search algorithm
#' obj_spec_random <- bekk_spec(init_values = "random")
#' x2 <- bekk_fit(obj_spec_random, StocksBonds, QML_t_ratios = FALSE,
#' max_iter = 50, crit = 1e-9)
#' summary(x1)
#' summary(x2)
#' plot(x1)
#' plot(x2)
#' # Fitting an asymmetric BEKK model with default starting values
#' obj_spec_fix <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE),
#' init_values = NULL)
#' x1 <- bekk_fit(obj_spec_fix, StocksBonds)
#' obj_spec_random <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE),
#' init_values = "random")
#' x2 <- bekk_fit(obj_spec_random, StocksBonds)
#' summary(x1)
#' summary(x2)
#' }
#' @export
bekk_spec <- function(model = list(type = "bekk", asymmetric = FALSE),
                      init_values = NULL, signs = NULL, N = NULL) {

  if(!is.logical(model$asymmetric) || is.null(model$asymmetric)){
    stop('Please specify whether the model to be estimated is asymmetric or not.')

  }
  # Checking type
  match.arg(model$type, c("bekk", "dbekk", "sbekk"))
  # Checking inputs
  if(!is.null(N) & is.numeric(init_values) ) {
    if(nrow(init_values) != 2 * N^2 + N * (N + 1)/2 & model$type == "bekk"  & model$asymmetric == FALSE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(nrow(init_values) != 3 * N^2 + N * (N + 1)/2 & model$type == "bekk"  & model$asymmetric == TRUE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(nrow(init_values) != 2 * N + N * (N + 1)/2 & model$type == "dbekk"  & model$asymmetric == FALSE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(nrow(init_values) != 3 * N + N * (N + 1)/2 & model$type == "dbekk"  & model$asymmetric == TRUE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(nrow(init_values) != 2 + N * (N + 1)/2 & model$type == "sbekk"  & model$asymmetric == FALSE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(nrow(init_values) != 3 + N * (N + 1)/2 & model$type == "sbekk"  & model$asymmetric == TRUE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
  }

  specification <- list(model = model, init_values = init_values, signs = signs, N = N)
  class(specification) <- 'bekkSpec'
  class(specification)[2] <- model$type

  if (model$asymmetric == TRUE) {

    if (!is.null(signs)) {
      signs <- matrix(signs, ncol = 1)
      if(any(abs(signs) != 1)){
        stop('Elements of "signs" must be either "1" or "-1".')
      }
    }

    specification$model$signs <- signs
    class(specification)[2] <- paste(class(specification)[2], 'a', sep = "")
  }

  return(specification)
}
