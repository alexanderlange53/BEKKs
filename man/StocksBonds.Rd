\name{StocksBonds}
\alias{StocksBonds}
\docType{data}
\title{
Daily stock and Bond returns
}
\description{
Bivariate data set consisting of daily S&P 500 bond and MSCI World returns from December 1995 to December 2019.
}
\usage{data("StocksBonds")}
\format{
  A data frame with 6073 observations on the following 2 variables.
  \describe{
    \item{S&P 500 Bonds}{a numeric vector}
    \item{MSCI World}{a numeric vector}
  }
}
\source{
Yahoo Finance.
}
\examples{
data(StocksBonds)
## maybe str(StocksBonds) ; plot(StocksBonds) ...
}
\keyword{datasets}
