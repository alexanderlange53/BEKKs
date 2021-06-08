BEKKs
=====

Modelling conditional volatilities and correlations among multivariate time series in R

## Overview

The package 'BEKKs' contains functions and methods for a careful analysis, estimation and forecasting of financial asset returns dynamics, and the construction and evaluation of financial portfolios. Modelling correlations and covarainces is important to determine portfolios with focus on hedging and asset specialization strategies, as well as to forecast value-at-risk (VaR) thresholds. 

The cornerstone functions to estimate correlation and covaraince processes are

-   `bekk()` the BEKK(1,1) model by Engle and Kroner (1995).

## Installation

```r
install.packages("BEKKs")
```

Alternatively, install the development version


```r
install.packages("devtools")
devtools::install_github("alexanderlange53/BEKKs")
```


```r
library("BEKKs")
```

## Usage

To get started, use the example data set which is included in the package. The data set consists of two daily financial time series, i.e. S&P 500 bond returns and MSCI World returns. More details on the data set are provided in the description file `?StocksBonds`.

```r
library("ggplot2")
library("ggfortify")
autoplot(StocksBonds  , facet = TRUE) + theme_bw()
```

![](man/figures/Data.png)

We estimate the conditional variance and covariance processes via the BEKK(1,1) model

```r
m1 <- bekk(StocksBonds)
```
