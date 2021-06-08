#' @export

summary.bekk <- function(object, ...) {
  bekkObject <- object
  cat(paste("\n", "BEKK estiamtion results", "\n", sep = ""))
  underScore <- paste(rep("-", nchar("BEKK estiamtion results")), collapse = "")
  cat(underScore)
  cat("\nLog-likelihood: ")
  cat(bekkObject$log_likelihood)
  cat("\nBEKK model stationary: ")
  cat(bekkObject$BEKK_valid)
  cat("\nNumber of BHHH iterations: ")
  cat(bekkObject$iter)
  cat("\nEstimated paramater matrices: \n")
  cat("\nC \n")
  print(bekkObject$C0)
  cat("\nA \n")
  print(bekkObject$A)
  cat("\nG \n")
  print(bekkObject$G)
  cat("\nt-values of paramater matrices: \n")
  cat("\nC \n")
  print(bekkObject$C0_t)
  cat("\nA \n")
  print(bekkObject$A_t)
  cat("\nG \n")
  print(bekkObject$G_t)
}
