VaR_BEKK<-function(theta,r,portfolio_weights,level){
  #compute H
  N=ncol(r)
  n=nrow(r)
  c_mats=coef_mat(theta, N)
  H=sigma_bekk(r, param_mat$c0, param_mat$a, param_mat$g)
  H_n1=c_mats$C0%*%t(c_mats$C0) + t(c_mats$a)t(r[n,])%*%r[n,]%*%c_mats$a+t(c_mats$g)%*%H[[n]]%*%c_mats$g
  nu=portfolio_weights%*%eigen_value_decomposition(H_n1)

  return(qnorm((1-level),sd=sum(portfolio_weights)))
}
