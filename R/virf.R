#' Estimating multivariate volatility impulse response functions (VIRF) for BEKK models
#'
#' @description Method for estimating VIRFs of N-dimensional BEKK models.
#'
#' @param fit An object of class "bekkfit" from function \link{bekk_fit}.
#' @param time Time instace to calculate VIRFs for.
#' @param q A vector specifying the quantiles to be considered for a shock on which basis the VIRFs are generated.
#' @param periods An integer defining the number periods for which the VIRFs are generated.
#' @return  Returns an object of class "bekkVIRF".
#'
#' @import xts
#' @import stats
#' @export

virf <- function(fit ,time = 1, q = 0.05, index_series = 1, periods = 10){

  if (!inherits(fit, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for fit')
  }

  if(!( periods%%1==0) || periods < 1){
    stop('Please provide a posive integer for periods')
  }

  if(!(index_series%%1==0) || index_series < 1){
    stop('Please provide a posive integer for index_series')
  }
  if(index_series > ncol(fit$data)){
    stop('Total number of indices in the data is lower than index_series')
  }
  if((time%%1!=0 || time < 1) && !inherits(time, 'Date')){
    stop('Please provide a posive integer or a date object for time')
  }else if(!(time%%1!=0 || time < 1) && time > nrow(fit$data)){
    stop('Total number of observations is exeded by time')
  }  else if(inherits(time, 'Date') && !is.numeric(fit$data[time])){
    stop('Provided date object is not included in data')
  }


  UseMethod('virf')

}

#' @export
virf.bekk <- function(fit, time = 1, q = 0.05, index_series=1, periods = 10) {

  N <- ncol(fit$data)
  data <- fit$data
  H <- matrix(fit$H_t[time,],N,N)
  #get quantiles of returns
  residuals = fit$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
  for(i in 1: N){
    if(i==index_series){
      shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
    }else{
      shocks[i] = sapply(0.5,FUN=quantile,x=as.matrix(residuals[,i]))

    }
  }

  VIRF = virf_bekk(H, fit$A, fit$G, matrix(shocks,ncol=N, nrow = 1), periods)
  dupl <- duplication_mat(N)
  elim <- elimination_mat(N)

  # for (i in 1:nrow(VIRF)) {
  #   tm <- matrix((dupl%*%VIRF[i,]), N, N, byrow = T)
  #   tm2 <- sqrt(solve(diag(abs(diag(tm)))))%*%tm%*%sqrt(solve(diag(abs(diag(tm)))))
  #   diag(tm2) <- sqrt(abs(diag(tm)))%*%solve(diag(abs(diag(tm))))%*%diag(diag(tm))
  #   VIRF[i,] <- elim%*%c(tm2)
  # }




  VIRF <- as.data.frame(VIRF)
  for(i in 1:ncol(VIRF)){
    colnames(VIRF)[i] <- paste("VIRF for", colnames(fit$sigma_t)[i], sep=" ")
    sub("Conditional","conditional",  colnames(VIRF)[i])
  }


  result <- list(VIRF=VIRF,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 fit=fit)
  class(result) <- c('virf','bekkFit', 'bekk')
  return(result)
}

virf.bekka <- function(fit, time = 1, q = 0.05, index_series=1, periods = 10) {

  N <- ncol(fit$data)
  H <- matrix(fit$H_t[time,],N,N)
  e
  #get quantiles of returns
  residuals = fit$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  for(i in 1: N){
    if(i==index_series){
    shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
    }else{
      shocks[i] = sapply(0.5,FUN=quantile,x=as.matrix(residuals[,i]))

    }
  }



  VIRF =  virf_bekka(H_t, fit$A, fit$B, fit$G, fit$signs, fit$expected_signs, shocks, periods)


  VIRF <- as.data.frame(VIRF)
  for(i in 1:ncol(VIRF)){
    colnames(VIRF)[i] <- paste("VIRF for", colnames(fit$sigma_t)[i], sep=" ")
    sub("Conditional","conditional",  colnames(VIRF)[i])
  }



  result <- list(VIRF=VIRF,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 fit=fit)
  class(result) <- c('virf','bekkFit', 'bekka')
  return(result)
}
