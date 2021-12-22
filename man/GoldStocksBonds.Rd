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
  The format is:
?zoo? series from 686361600 to 1632873600
  Data: num [1:7346, 1:3] 0.00578 0.00224 -0.0014 0.0095 -0.00586 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:3] "Gold" "S&P 500" "US Treasury Bond Future"
  Index:  num [1:7346] 6.86e+08 6.86e+08 6.87e+08 6.87e+08 6.87e+08 ...
 - attr(*, "tzone")= chr "UTC"
 - attr(*, "tclass")= chr "Date"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(GoldStocksBonds)
## maybe str(GoldStocksBonds) ; plot(GoldStocksBonds) ...
}
\keyword{datasets}
