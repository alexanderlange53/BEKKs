#' Estimating multivariate volatility impulse response functions (VIRF) for BEKK models
#'
#' @description Method for estimating VIRFs of N-dimensional BEKK models.
#'
#' @param x An object of class "bekkfit" from function \link{bekk_fit}.
#' @param time Time instace to calculate VIRFs for.
#' @param q A vector specifying the quantiles to be considered for a shock on which basis the VIRFs are generated.
#' @param n.ahead An integer defining the number periods for which the VIRFs are generated.
#' @return  Returns an object of class "bekkVIRF".
#'
#' @import xts
#' @import stats
#' @export

virf <- function(x ,time = 1, q = 0.05, index_series = 1, n.ahead = 10, ci = 0.9){

  if (!inherits(x, 'bekkFit')) {
    stop('Please provide and object of class "bekkFit" for x')
  }

  if(!( n.ahead%%1==0) || n.ahead < 1){
    stop('Please provide a posive integer for periods')
  }

  if(!(index_series%%1==0) || index_series < 1){
    stop('Please provide a posive integer for index_series')
  }
  if(index_series > ncol(x$data)){
    stop('Total number of indices in the data is lower than index_series')
  }
  if((time%%1!=0 || time < 1) && !inherits(time, 'Date')){
    stop('Please provide a posive integer or a date object for time')
  }else if(!(time%%1!=0 || time < 1) && time > nrow(x$data)){
    stop('Total number of observations is exeded by time')
  }  else if(inherits(time, 'Date') && !is.numeric(x$data[time])){
    stop('Provided date object is not included in data')
  }


  UseMethod('virf')

}

#' @export
virf.bekk <- function(x, time = 1, q = 0.05, index_series=1, n.ahead = 10, ci = 0.9) {

  N <- ncol(x$data)
  data <- x$data
  H <- matrix(x$H_t[time,],N,N)
  #get quantiles of returns
  residuals = x$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
  for(i in 1: N){
    if(i==index_series){
      shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
    }else{
      shocks[i] = sapply(0.5,FUN=quantile,x=as.matrix(residuals[,i])) * 0

    }
  }

  VIRF = virf_bekk(H, x$theta, matrix(shocks, ncol=N, nrow = 1), n.ahead)
  #dupl <- duplication_mat(N)
  #elim <- elimination_mat(N)

  score_final = score_bekk(x$theta, x$data)
  s1_temp = solve(t(score_final) %*% score_final)
  s1 = eigen_value_decomposition(s1_temp)
  print(s1)
  theta = x$theta
  theta_lower = theta - s1%*%as.matrix(rep(qnorm(ci), nrow(x$theta)))
  x_lower <- coef_mat(theta_lower, N)
  H_lower = matrix(sigma_bekk(x$data,x_lower$c0, x_lower$a, x_lower$g)$sigma_t[time,], N, N)
  VIRF_lower = virf_bekk(H_lower, theta_lower, matrix(shocks,ncol=N, nrow = 1), n.ahead)



  theta_upper = theta + s1%*%as.matrix(rep(qnorm(ci), nrow(x$theta)))
  x_upper <- coef_mat(theta_upper, N)
  H_upper = matrix(sigma_bekk(x$data, x_upper$c0, x_upper$a, x_upper$g)$sigma_t[time,], N, N)

  VIRF_upper = virf_bekk(H_upper, theta_upper, matrix(shocks,ncol=N, nrow = 1), n.ahead)

  # for (i in 1:nrow(VIRF)) {
  #   tm <- matrix((dupl%*%VIRF[i,]), N, N, byrow = T)
  #   tm2 <- sqrt(solve(diag(abs(diag(tm)))))%*%tm%*%sqrt(solve(diag(abs(diag(tm)))))
  #   diag(tm2) <- sqrt(abs(diag(tm)))%*%solve(diag(abs(diag(tm))))%*%diag(diag(tm))
  #   VIRF[i,] <- elim%*%c(tm2)
  # }




  VIRF <- as.data.frame(VIRF)
  VIRF_lower <- as.data.frame(VIRF_lower)
  VIRF_upper <- as.data.frame(VIRF_upper)
  for(i in 1:ncol(VIRF)){
    colnames(VIRF)[i] <- paste("VIRF for", colnames(x$sigma_t)[i], sep=" ")
    sub("Conditional","conditional",  colnames(VIRF)[i])
  }



  result <- list(VIRF=VIRF,
                 VIRF_upper=VIRF_upper,
                 VIRF_lower=VIRF_lower,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 x=x)
  class(result) <- c('virf','bekkFit', 'bekk')
  return(result)
}

virf.bekka <- function(x, time = 1, q = 0.05, index_series=1, n.ahead = 10) {

  N <- ncol(x$data)
  H <- matrix(x$H_t[time,],N,N)
  e
  #get quantiles of returns
  residuals = x$e_t
  shocks = matrix(0, nrow = 1, ncol = N)

  for(i in 1: N){
    if(i==index_series){
    shocks[index_series] = sapply(q,FUN=quantile,x=as.matrix(residuals[,index_series]))
    }else{
      shocks[i] = sapply(0.5,FUN=quantile,x=as.matrix(residuals[,i]))

    }
  }



  VIRF =  virf_bekka(H_t, t(x$C0), x$A, x$B, x$G, x$signs, x$expected_signs, shocks, n.ahead)


  VIRF <- as.data.frame(VIRF)
  for(i in 1:ncol(VIRF)){
    colnames(VIRF)[i] <- paste("VIRF for", colnames(x$sigma_t)[i], sep=" ")
    sub("Conditional","conditional",  colnames(VIRF)[i])
  }



  result <- list(VIRF=VIRF,
                 N=N,
                 time=time,
                 q=q,
                 index_series=index_series,
                 x=x)
  class(result) <- c('virf','bekkFit', 'bekka')
  return(result)
}
