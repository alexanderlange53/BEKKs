#' Calculating Value-at-Risk (VaR)
#'
#' @description Method for calculating VaR from estimated covariance processes (\link{bekk_fit}) or predicted covariances (\link{predict}).
#'
#' @param x An object of class "bekkFit" from the function \link{bekk_fit} or an object of class "bekkForecast" from the function \link{predict}.
#' @param p A numerical value that determines the confidence level. The default value is set at 0.99 in accordance with the Basel Regulation.
#' @param portfolio_weights A vector determining the portfolio weights to calculate the portfolio VaR. If set to "NULL", the univariate VaR for each series are calculated.
#' @param distribution A character string determining the assumed distribution of the residuals. Implemented are "normal", "empirical" and "t". The default is using the empirical distribution of the residuals.
#' @return  Returns a S3 class "var" object containing the VaR forecast and respective confidence bands.
#' @examples
#' \donttest{
#'
#' data(StocksBonds)
#' obj_spec <- bekk_spec()
#' x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' # single VaRs of series
#' x2 <- VaR(x1, distribution="normal")
#' plot(x2)
#'
#' # VaR of equally-weighted portfolio
#' portfolio_weights <- c(0.5, 0.5)
#' x3 <- VaR(x1, portfolio_weights = portfolio_weights)
#' plot(x3)
#'
#' # VaR of traditional 30/70 weighted bond and stock portfolio
#' portfolio_weights <- c(0.3, 0.7)
#' x4 <- VaR(x1, portfolio_weights = portfolio_weights)
#' plot(x4)
#'
#' }
#'
#' @import xts
#' @import stats
#' @import moments
#' @export

VaR <- function(x, p = 0.99, portfolio_weights = NULL, distribution = "empirical") {

  UseMethod('VaR')
}

#' @export
VaR.bekkFit <-  function(x, p = 0.99, portfolio_weights = NULL, distribution = "empirical" )
{
  if(nrow(x$data) < 1000 && distribution == "empirical"){
    stop("Using the empirical distribution is not stable for time series with less than 1000 observations!")
  }

  alpha = p
  N = ncol(x$data)
  n = nrow(x$data)
  match.arg(distribution, c("empirical", "t", "normal"))

  if(length(portfolio_weights)!= N && !is.null(portfolio_weights)){
    stop("Portfolio weights do no match number of time series")
  }


  #specify quantiles here
  if(distribution == "t"){
    #fit skewed t
    skew_t <- function(i){
      kurtos = moments::kurtosis(x$e_t[,i])-3
      df = 6/kurtos+4
      if(df <= 4){
        df=4.001
      }
      return(qt(1-alpha, df = df)/sqrt(df/(df-2)))
    }

    qtls <- sapply(1:N, skew_t)
  }else if(distribution == "empirical"){
    empirical <- function(i){
      quantile(x$e_t[,i],1-alpha)
    }
    qtls <- sapply(1:N, empirical)

  } else if(distribution == "normal"){
    #fit skewed t

    qtls <- rep(qnorm(1-alpha),ncol(x$data))
  } else{
    qtls <- rep(qnorm(1-alpha),ncol(x$data))
  }
  #quantile(x$e_t, )

  if (is.null(portfolio_weights)) {
    columns = ncol(x$data)
    csd <- extract_csd(x)
    VaR <- matrix(NA, nrow = nrow(x$data), ncol = ncol(x$data))

    for(i in 1:n) {
       for(j in 1: ncol(x$data)){
      VaR[i, j] =  sqrt(matrix(x$H_t[i,],N,N)[j,j]) * qtls[j]
        }
    }
    VaR <- as.data.frame(VaR)
    for(column in 1:columns) {
      r = as.vector(na.omit(x$data[,column]))
      if (!is.numeric(r)) stop("The selected column is not numeric")
      m2 =  csd[, column]
      #VaR[, column] = - qnorm(alpha)*m2
      #VaR <- as.data.frame(VaR)

      for (i in 1:ncol(x$data)) {
        colnames(VaR)[i] <- paste('VaR of', colnames(x$data)[i])
      }
    }
  } else {
    if(distribution == "t"){
      #fit skewed t

        kurtos = moments::kurtosis(x$e_t%*%portfolio_weights)-3
        df = 6/kurtos+4
        if(df <= 4){
          df=4.001
        }
        qtls <- qt(1-alpha, df = df)/sqrt(df/(df-2))
      }
    else if(distribution == "empirical"){

      qtls <- quantile(x$e_t%*%portfolio_weights,1-alpha)

      } else if(distribution == "normal"){
      #fit skewed t

      qtls <-qnorm(1-alpha)
    } else{
      qtls <- qnorm(1-alpha)
    }


    VaR <- matrix(NA, nrow = nrow(x$data), ncol = 1)
    for(i in 1:nrow(x$H_t)) {
      VaR[i,] <- qtls*sqrt(portfolio_weights%*%matrix(x$H_t[i,], ncol = ncol(x$data))%*%portfolio_weights)
      #VaR[i,] <- portfolio_weights%*%eigen_value_decomposition(matrix(x$H_t[i,], ncol = ncol(x$data)))%*%qtls
      }
    VaR <- as.data.frame(VaR)
  }

  if (inherits(x$data, "ts")) {
    VaR <- ts(VaR, start = time(x$data)[1], frequency = frequency(x$data))
  }else if(inherits(x$data, "xts") || inherits(x$data, "zoo") ){
    VaR <- xts(VaR, order.by = time(x$data[1:nrow(x$data),]))
      }

  result <- list(VaR = VaR,
                 p = p,
                 portfolio_weights = portfolio_weights,
                 bekk = x)
  class(result) <- c('var', 'bekkFit')
  return(result)
}


#' @export
VaR.bekkForecast <-  function(x, p = 0.99, portfolio_weights = NULL, distribution = "empirical")
{

  if(nrow(x$bekkfit$data) < 1000 && distribution == "empirical"){
    stop("Using the empirical distribution is not stable for time series with less than 1000 observations!")
  }

  alpha = p


  obj <- x$bekkfit
  obj$H_t <- rbind(x$bekkfit$H_t, x$H_t_forecast)
  obj$sigma_t <- rbind(x$bekkfit$sigma_t, x$volatility_forecast)
  N = ncol(obj$data)
  n = nrow(obj$H_t)
  if(length(portfolio_weights)!= N && !is.null(portfolio_weights)){
    stop("Portfolio weights do no match number of time series")
  }

  if(distribution == "t"){
    #fit skewed t
    skew_t <- function(i){
      kurtos = moments::kurtosis(obj$e_t[,i])-3
      df = 6/kurtos+4
      if(df <= 4){
        df=4.001
      }
      return(qt(1-alpha, df = df)/sqrt(df/(df-2)))
    }

    qtls <- sapply(1:N, skew_t)
  }else if(distribution == "empirical"){
    empirical <- function(i){
      quantile(obj$e_t[,i],1-alpha)
    }
    qtls <- sapply(1:N, empirical)

  } else if(distribution == "normal"){
    #fit skewed t

    qtls <- rep(qnorm(1-alpha),ncol(obj$data))
  } else{
    qtls <- rep(qnorm(1-alpha),ncol(obj$data))
  }

  if (is.null(portfolio_weights)) {
    columns = ncol(x$bekkfit$data)
    #csd <- extract_csd(obj)
    VaR <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = ncol(x$bekkfit$data))


    for(i in 1:nrow(obj$H_t)) {
      for(j in 1: ncol(x$bekkfit$data)){
        VaR[i, j] =  sqrt(matrix(obj$H_t[i,],N,N)[j,j]) * qtls[j]
      }
    }
    VaR <- as.data.frame(VaR)

      for (i in 1:ncol(x$bekkfit$data)) {
        colnames(VaR)[i] <- paste('VaR of', colnames(x$bekkfit$data)[i])
      }


    # Confidence intervals
    H_t_lower <- rbind(x$bekkfit$H_t[-nrow(x$bekkfit$H_t),], x$H_t_lower_conf_band)
    H_t_upper <- rbind(x$bekkfit$H_t[-nrow(x$bekkfit$H_t),], x$H_t_upper_conf_band)

    colnames(x$volatility_lower_conf_band) = colnames(x$volatility_forecast)
    obj$sigma_t <- rbind(x$bekkfit$sigma_t, x$volatility_lower_conf_band)

    csd_lower <- extract_csd(obj)
    VaR_lower <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = ncol(x$bekkfit$data))

    for(i in 1:nrow(obj$H_t)) {
      for(j in 1: ncol(x$bekkfit$data)){
        VaR_lower[i, j] =  sqrt(matrix(H_t_lower[i,],N,N)[j,j]) * qtls[j]
      }
    }

      VaR_lower <- as.data.frame(VaR_lower)

      for (i in 1:ncol(x$bekkfit$data)) {
        colnames(VaR_lower)[i] <- paste('VaR of', colnames(x$bekkfit$data)[i])
      }

    colnames(x$volatility_upper_conf_band) = colnames(x$volatility_forecast)
    obj$sigma_t <- rbind(x$bekkfit$sigma_t, x$volatility_upper_conf_band)

    csd_upper <- extract_csd(obj)
    VaR_upper <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = ncol(x$bekkfit$data))
    for(i in 1:nrow(obj$H_t)) {
      for(j in 1: ncol(x$bekkfit$data)){
        VaR_upper[i, j] =  sqrt(matrix(H_t_upper[i,],N,N)[j,j]) * qtls[j]
      }
    }

      VaR_upper <- as.data.frame(VaR_upper)

      for (i in 1:ncol(x$bekkfit$data)) {
        colnames(VaR_upper)[i] <- paste('VaR of', colnames(x$bekkfit$data)[i])
      }



  } else {

    if(distribution == "t"){
      #fit skewed t

      kurtos = moments::kurtosis(obj$e_t%*%portfolio_weights)-3
      df = 6/kurtos+4
      if(df <= 4){
        df=4.001
      }
      qtls <- qt(1-alpha, df = df)/sqrt(df/(df-2))
    }
    else if(distribution == "empirical"){

      qtls <- quantile(obj$e_t%*%portfolio_weights,1-alpha)

    } else if(distribution == "normal"){
      #fit skewed t

      qtls <-qnorm(1-alpha)
    } else{
      qtls <- qnorm(1-alpha)
    }

    VaR <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = 1)

    for(i in 1:nrow(obj$H_t)) {
      #VaR[i,] <- -qnorm(alpha)*sqrt(portfolio_weights%*%matrix(x$H_t[i,], ncol = ncol(x$data))%*%portfolio_weights)
      VaR[i,] <- sqrt(portfolio_weights%*%matrix(obj$H_t[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)*qtls
    }
    # for(i in 1:nrow(obj$H_t)) {
    #   VaR[i,] <- -qnorm(alpha)*sqrt(portfolio_weights%*%matrix(obj$H_t[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)
    # }
    VaR <- as.data.frame(VaR)

    # Confidence intervals
    VaR_lower <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = 1)
    VaR_upper <- matrix(NA, nrow = nrow(x$bekkfit$data) + x$n.ahead, ncol = 1)

    H_t_lower <- rbind(x$bekkfit$H_t[-nrow(x$bekkfit$H_t),], x$H_t_lower_conf_band)
    H_t_upper <- rbind(x$bekkfit$H_t[-nrow(x$bekkfit$H_t),], x$H_t_upper_conf_band)

    # for(i in 1:nrow(obj$H_t)) {
    #   VaR_lower[i,] <- -qnorm(alpha)*sqrt(portfolio_weights%*%matrix(H_t_lower[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)
    #   VaR_upper[i,] <- -qnorm(alpha)*sqrt(portfolio_weights%*%matrix(H_t_upper[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)
    # }

    for(i in 1:nrow(obj$H_t)) {
      VaR_lower[i,] <- sqrt(portfolio_weights%*%matrix(H_t_lower[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)*qtls
      VaR_upper[i,] <- sqrt(portfolio_weights%*%matrix(H_t_upper[i,], ncol = ncol(x$bekkfit$data))%*%portfolio_weights)*qtls
    }
    VaR_lower <- as.data.frame(VaR_lower)
    VaR_upper <- as.data.frame(VaR_upper)
  }

  if (inherits(x$data, "ts")) {
    VaR <- ts(VaR, start = time(x$bekkfit$data)[1], frequency = frequency(x$bekkfit$data))
  }  else if(inherits(x$data, "xts") || inherits(x$data, "zoo") ){
    VaR <- xts(VaR, order.by = time(x$data[1:nrow(x$data),]))
  }

  result <- list(VaR = VaR,
                 VaR_lower = VaR_lower,
                 VaR_upper = VaR_upper,
                 p = p,
                 portfolio_weights = portfolio_weights,
                 n.ahead = x$n.ahead,
                 bekk = x)
  class(result) <- c('var', 'bekkForecast')
  return(result)
}
