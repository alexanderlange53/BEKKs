\name{GoldStocksBonds}
\alias{GoldStocksBonds}
\docType{data}
\title{
Gold stock and Bond returns
}
\description{
Trivariate data set consisting of daily gold, S&P 500 and U.S. Treasury Bond Future returns from October 1991 to October 2021.
}
\usage{data("GoldStocksBonds")}
\format{
  A time series matrix of class mts 7346 observations on the following 3 variables.
  \describe{
    \item{Gold}{a numeric vector}
    \item{S&P 500}{a numeric vector}
    \item{US Treasury Bond Future}{a numeric vector}
  }
}

\source{
Yahoo Finance.
}
\examples{
data(GoldStocksBonds)
## maybe str(GoldStocksBonds) ; plot(GoldStocksBonds) ...
}
\keyword{datasets}
