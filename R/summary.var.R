#' @export
print.var <- function(x,...){
  object <- x
  bekkObject <- object$bekk
  if (any(class(bekkObject) == 'bekk')) {
    cat(paste("\n", "BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("BEKK VaR results")), collapse = "")
  } else if (any(class(bekkObject) == 'bekka')) {
    cat(paste("\n", "Asymmetric BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric BEKK VaR results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekk')) {
    cat(paste("\n", "Diagonal BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Diagonal BEKK VaR results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekka')) {
    cat(paste("\n", "Asymmetric diagonal BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric diagonal BEKK VaR results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekk')) {
    cat(paste("\n", "Scalar BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Scalar BEKK VaR results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekka')) {
    cat(paste("\n", "Asymmetric scalar BEKK VaR results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric scalar BEKK VaR results")), collapse = "")
  }
  cat(underScore)
  cat("\nValue-at-risk confidence level: ")
  cat(as.character(object$p))
  if(!is.null(object$portfolio_weights)){
    cat("\nPortfolio weights: ")
    cat(object$portfolio_weights)
  }else{
    cat("\nPortfolio weights: None\n")

  }
}
