#' @export

summary.backtest <- function(object, ...) {
  bekkObject <- object$bekkFit
  if (any(class(bekkObject) == 'bekk')) {
    cat(paste("\n", "BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'bekka')) {
    cat(paste("\n", "Asymmetric BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekk')) {
    cat(paste("\n", "Diagonal BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Diagonal BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekka')) {
    cat(paste("\n", "Asymmetric diagonal BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric diagonal BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekk')) {
    cat(paste("\n", "Scalar BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Scalar BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekka')) {
    cat(paste("\n", "Asymmetric scalar BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric scalar BEKK backtesting results")), collapse = "")
  }
  cat(underScore)
  cat("\nValue-at-risk confidence level: ")
  cat(as.character(object$p))
  cat("\nWindow length: ")
  cat(as.character(object$window_length))
  if(!is.null(object$portfolio_weights)){
    res_hit <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=1))
    res_Kupiec <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=2))
    res_Christoffesen <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=2))
    colnames(res_hit)=c("")
    colnames(res_Kupiec)=c("")
    colnames(res_Christoffesen)=c("")
    row.names(res_hit)=c("")
    row.names(res_Kupiec)=c("Test", "p-value")
    row.names(res_Christoffesen)=c("Test", "p-value")

      res_Kupiec[,1]=object$backtests$LRuc
      res_Christoffesen[,1]=object$backtests$LRcc


    cat("\nPortfolio weights: ")
    cat(object$portfolio_weights)
    cat("\n")
    cat(underScore)
    res_hit[1,]=object$hit_rate
    cat(paste("\nHit rate:", round(res_hit,3), "\n", sep = " " ))

    cat("\nUnconditional coverage test of Kupiec: \n")
    print(res_Kupiec)
    cat("\nconditional coverage test of Christoffesen: \n")
    print(res_Christoffesen)
  }else{
    res_hit <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=1))
    res_Kupiec <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=2))
    res_Christoffesen <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=2))
    colnames(res_hit)=colnames(object$out_sample_returns)
    colnames(res_Kupiec)=colnames(res_hit)
    colnames(res_Christoffesen)=colnames(res_hit)
    row.names(res_hit)=c("")
    row.names(res_Kupiec)=c("Test", "p-value")
    row.names(res_Christoffesen)=c("Test", "p-value")
    for(i in 1:ncol(object$VaR)){
      res_Kupiec[,i]=object$backtests[[i]]$LRuc
      res_Christoffesen[,i]=object$backtests[[i]]$LRcc
    }

    cat("\nPortfolio weights: None\n")
    cat(underScore)
    res_hit[1,]=object$hit_rate
    cat("\nHit rates: \n")
    cat(underScore)
    print(res_hit)
    cat("\nUnconditional coverage test of Kupiec: \n")
    print(res_Kupiec)
    cat("\nconditional coverage test of Christoffesen: \n")
    print(res_Christoffesen)
  }

}

#' @export
print.backtest <- function(x,...){
  object <- x
  bekkObject <- object$bekkFit
  if (any(class(bekkObject) == 'bekk')) {
    cat(paste("\n", "BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'bekka')) {
    cat(paste("\n", "Asymmetric BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekk')) {
    cat(paste("\n", "Diagonal BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Diagonal BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'dbekka')) {
    cat(paste("\n", "Asymmetric diagonal BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric diagonal BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekk')) {
    cat(paste("\n", "Scalar BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Scalar BEKK backtesting results")), collapse = "")
  } else if (any(class(bekkObject) == 'sbekka')) {
    cat(paste("\n", "Asymmetric scalar BEKK backtesting results", "\n", sep = ""))
    underScore <- paste(rep("-", nchar("Asymmetric scalar BEKK backtesting results")), collapse = "")
  }
  cat(underScore)
  cat("\nValue-at-risk confidence level: ")
  cat(as.character(object$p))
  cat("\nWindow length: ")
  cat(as.character(object$window_length))
  if(!is.null(object$portfolio_weights)){
    cat("\nPortfolio weights: ")
    cat(object$portfolio_weights)
    cat("\n")
    cat(underScore)
    cat("\nHit rate: ")
    cat(round(object$hit_rate,3))
  }else{
    cat("\nPortfolio weights: None\n")
    res_hit <- data.frame(matrix(NA, ncol = ncol(object$VaR), nrow=1))
    colnames(res_hit)=colnames(object$out_sample_returns)
    row.names(res_hit)=c("")

    res_hit[1,]=round(object$hit_rate,3)
    cat(underScore)
    cat("\nHit rates: \n")
    print(res_hit)

  }
}
