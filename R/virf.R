#' Estimating multivariate volatility impulse response functions (VIRF) for BEKK models
#'
#' @description Method for estimating VIRFs of N-dimensional BEKK models. Currently, only VIRFs for symmetric BEKK models are implemented.
#'
#' @param x An object of class "bekkfit" from function \link{bekk_fit}.
#' @param time Time instace to calculate VIRFs for.
#' @param q A vector specifying the quantiles to be considered for a shock on which basis the VIRFs are generated.
#' @param n.ahead An integer defining the number periods for which the VIRFs are generated.
#' @return  Returns an object of class "bekkVIRF".
#'
#' @import xts
#' @import stats
#' @import numDeriv
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

  # score_final = score_bekk(x$theta, x$data)
  # s1_temp = solve(t(score_final) %*% score_final)
  # s1 = eigen_value_decomposition(s1_temp)
  hesse_final = solve(hesse_bekk(x$theta, x$data))
  #s1_temp = solve(hesse_final)

  s1_temp = function(th){
    virf_bekk(H, th, matrix(shocks, ncol=N, nrow = 1), n.ahead)
  }

  th<-x$theta
  d_virf = jacobian(s1_temp,th)
  s1_temp=d_virf%*%hesse_final%*%t(d_virf)

  print(s1_temp)
#   s1 = s1_temp*0
#   counter = 1
#   while(counter < nrow(s1)){
#     s1[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]=s1_temp[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]
#   counter = counter + n.ahead
# }

  s1 = sqrt(diag(s1_temp)) * qnorm(ci)
  #return(s1)
  #print(det(d_virf%*%hesse_final%*%t(d_virf)))
  VIRF_lower = VIRF  - matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

  VIRF_upper = VIRF + matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

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
#' @export
virf.dbekk <- function(x, time = 1, q = 0.05, index_series=1, n.ahead = 10, ci = 0.9) {

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

  VIRF = virf_dbekk(H, x$theta, matrix(shocks, ncol=N, nrow = 1), n.ahead)

  hesse_final = solve(hesse_dbekk(x$theta, x$data))

  s1_temp = function(th){
    virf_bekk(H, th, matrix(shocks, ncol=N, nrow = 1), n.ahead)
  }

  th<-x$theta
  d_virf = jacobian(s1_temp,th)
  s1_temp=d_virf%*%hesse_final%*%t(d_virf)
  #   s1 = s1_temp*0
  #   counter = 1
  #   while(counter < nrow(s1)){
  #     s1[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]=s1_temp[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]
  #   counter = counter + n.ahead
  # }

  s1 = sqrt(diag(s1_temp)) * qnorm(ci)
  #return(s1)
  #print(det(d_virf%*%hesse_final%*%t(d_virf)))
  VIRF_lower = VIRF  - matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

  VIRF_upper = VIRF + matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

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
#' @export
virf.sbekk <- function(x, time = 1, q = 0.05, index_series=1, n.ahead = 10, ci = 0.9) {

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

  VIRF = virf_sbekk(H, x$theta, matrix(shocks, ncol=N, nrow = 1), n.ahead)
  #dupl <- duplication_mat(N)
  #elim <- elimination_mat(N)

  # score_final = score_bekk(x$theta, x$data)
  # s1_temp = solve(t(score_final) %*% score_final)
  # s1 = eigen_value_decomposition(s1_temp)
  hesse_final = solve(hesse_scalar_bekk(x$theta, x$data))
  #s1_temp = solve(hesse_final)

  s1_temp = function(th){
    virf_sbekk(H, th, matrix(shocks, ncol=N, nrow = 1), n.ahead)
  }

  th<-x$theta
  d_virf = jacobian(s1_temp,th)
  s1_temp=d_virf%*%hesse_final%*%t(d_virf)
  #   s1 = s1_temp*0
  #   counter = 1
  #   while(counter < nrow(s1)){
  #     s1[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]=s1_temp[counter:(counter+n.ahead-1),counter:(counter+n.ahead-1)]
  #   counter = counter + n.ahead
  # }

  s1 = sqrt(diag(s1_temp)) * qnorm(ci)  #return(s1)
  #print(det(d_virf%*%hesse_final%*%t(d_virf)))
  VIRF_lower = VIRF  - matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

  VIRF_upper = VIRF + matrix(s1, nrow = n.ahead, ncol = N*(N+1)/2)

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
  class(result) <- c('virf','bekkFit', 'sbekk')
  return(result)
}
