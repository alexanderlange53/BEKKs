#' BEKK specification method
#'
#' @param model List containing the model type specification: Currently implemented are "bekk", "dbekk", "sbekk".
#' Moreover it can be specified whether the model should be estimated allowing for asymmetric volatility structure.
#' @param init_values initial values for \href{bekk_fit} during BHHH algorithm. It can be either a numerical vector of suitable dimension, or a character vector i.e. "random" to use a random starting values generator, or
#'  "simple" for relying on a simple initial values generator based on typical values for BEKK parameter found in the literature. If object from this function is passed to \link{bekk_sim}, init_values are used as parameters for data generating process.
#' @param signs Vector specifying asymmetry.
#' @param N Integer specifying the dimension of the BEKK model. Only relevant for \link{bekk_sim}.
#'
#'
#'
#'
#' @export


bekk_spec <- function(model = list(type = "bekk", asymmetric = FALSE),
                      init_values = NULL, signs = NULL, N = NULL, compare=FALSE) {


  # Checking inputs
  if(!is.null(N) & is.numeric(init_values) ) {
    if(length(init_values) != 2 * N^2 + N * (N + 1)/2 & type == "bekk"  & asymmetric = FALSE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
    if(length(init_values) != 3 * N^2 + N * (N + 1)/2 & type == "bekk"  & asymmetric = TRUE) {
      stop('Number of initial parameter does not match dimension of data and model.')
    }
  }

  specification <- list(model = model, init_values = init_values, signs = signs, N = N)
  class(specification) <- 'bekkSpec'
  class(specification)[2] <- model$type

  if (model$asymmetric == TRUE) {
    if(is.null(model$signs)){
      model$signs = rep(-1,N)
    }
    if(any(abs(signs)!=1)){
      stop('Elements of "signs" must be either "1" or "-1".')
    }
    if(length(model$signs)!=N){
      stop('Length of "signs" does not match dimension of data.')
    }
    class(specification)[2] <- paste(class(specification)[2], 'a', sep = "")
  }

  return(specification)

}
