VaR <- function(x) {
  UseMethod('VaR')
}

VaR.bekkFit =  function(x, p = 0.95)
{
  alpha = p#.setalphaprob(p)

  columns = ncol(x$data)
  csd <- extract_csd(x)
  VaR <- matrix(NA, nrow = nrow(x$data), ncol = ncol(x$data))

  for(column in 1:columns) {
    r = as.vector(na.omit(x$data[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric")
    m2 =  csd[, column] #centeredmoment(r,2)
    VaR[, column] = - mean(r) - qnorm(alpha)*m2
  }

  return(VaR)
}

VaR_BEKK <- function(theta, r, portfolio_weights){
  #compute H
  nc = ncol(r)
  n = nrow(r)

  c_mats = coef_mat(theta, nc)
  H = sigma_bekk(r, param_mat$c0, param_mat$a, param_mat$g)
  H_n1 = c_mats$C0%*%t(c_mats$C0) + t(c_mats$a)%*%t(r[n,])%*%r[n,]%*%c_mats$a+t(c_mats$g)%*%H[[n]]%*%c_mats$g
  nu = portfolio_weights%*%eigen_value_decomposition(H_n1)

  return(qnorm((1-nu), sd=sum(portfolio_weights)))
}
