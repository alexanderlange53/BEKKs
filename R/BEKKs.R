#' BEKKs: Volatility modelling
#'
#' @docType package
#' @name BEKKs
#' @author \itemize{
#' \item Markus Fülle  \email{fuelle@uni-goettingen.de}
#' \item Helmut Herwartz \email{hherwartz@uni-goettingen.de}
#' \item Alexander Lange \email{alexander.lange@uni-goettingen.de}
#' }
#' @description
#' This package implements MGARCH estimation techniques for conditional volatility modelling.\cr
#' @details
#' The main functions to retrieve structural impact matrices are:
#' \itemize{
#' \item \tabular{ll}{ \code{\link{bekk}} \tab Estimates a BEKK(1,1,1) model,}
#' }
#' @useDynLib BEKKs
#' @importFrom Rcpp sourceCpp
NULL
