# bhh_bekk <- function(r, theta, max_iter, crit) {
#
#   steps <- c(5,2,1,0.5,0.25,0.1,0.01,0.005,0)
#   step <- 0.01
#   count_loop <- 1
#   theta_candidate <- theta
#   exit_loop <- 0
#
#
#   while (count_loop < max_iter & exit_loop == 0) {
#     theta <- theta_candidate
#     theta_temp <- matrix(0, nrow = nrow(theta), length(steps))
#
#     score_function <- score_bekk(theta, r)
#     outer_score <- crossprod(score_function)
#     score_function <- colSums(score_function)
#
#     # Hier likelihood
#     lik <- loglike_bekk(theta, r)
#
#     for (i in 1:length(steps)) {
#       temp <- theta_candidate + step * steps[i]*solve(outer_score) %*% score_function # hier nochmal gucken H2 hat da kontrolliert on inverser exstiert
#       theta_temp[, i] <- temp
#     }
#
#     likelihood_candidates <- rep(0, length(steps))
#     likelihood_candidates[length(steps)] <- lik
#
#     j <- length(steps) - 1
#     exit_inner <-0
#     while (j >= 1 & exit_inner == 0) {
#       likelihood_candidates[j] <- loglike_bekk(theta_temp[, j], r)
#       if (likelihood_candidates[j+1] > likelihood_candidates[j]) {
#         exit_inner <-1
#       }
#       j <- j-1
#     }
#
#     max_index <- which.max(likelihood_candidates[(j+1):length(likelihood_candidates)]) + j
#     likelihood_best <- likelihood_candidates[max_index]
#
#     # exit criterion strange
#     if ((likelihood_best - likelihood_candidates[length(likelihood_candidates)])^2/abs(likelihood_candidates[length(likelihood_candidates)]) < crit) {
#       exit_loop <- 1
#     }
#     theta_candidate <- matrix(theta_temp[, max_index], ncol= 1)
#     count_loop <- count_loop + 1
#   }
#
#   likelihood_final <- loglike_bekk(theta_candidate, r)
#   score_final <- score_bekk(theta_candidate, r)
#   s1 <- sqrt(diag(solve(crossprod(score_final))))
#
#   t_val <- theta_candidate/s1
#
#   return(list(theta = theta_candidate,
#               t_val = t_val,
#               likelihood = likelihood_final,
#               iter = count_loop))
# }
#
# # score_bekk <- function(theta, r) {
# #   N <- ncol(r)
# #   N2 <- N^2
# #
# #   L_elimination <- elimination_mat(N)
# #   D_duplication <- duplication_mat(N)
# #
# #   D_gen_inv <- solve(crossprod(D_duplication)) %*% t(D_duplication)
# #
# #   gradients <- matrix(0, nrow = nrow(r), ncol = nrow(theta))
# #
# #   c_mat <-  coef_mat(theta, N)
# #
# #   c0 <- c_mat$c0
# #   a <- c_mat$a
# #   g <- c_mat$g
# #
# #   c_full <- crossprod(c0)
# #
# #   # Partial derivatives for initial period t = 1
# #   ht <- crossprod(r)/nrow(r)
# #
# #   dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% ht)
# #   dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% ht)
# #   dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination)
# #
# #   dHdtheta <- t(cbind(dHdc, dHda, dHdg))
# #
# #   ht_sqrt_inv <- solve(sqrtm(ht))
# #
# #   et <- ht_sqrt_inv %*% r[1,]
# #   # Score function
# #   for (k in 1:nrow(theta)) {
# #     dh <- matrix(dHdtheta[k, ], N, N)
# #
# #     mat_temp <- ht_sqrt_inv %*% dh %*% ht_sqrt_inv %*% (diag(N) - tcrossprod(et))
# #
# #     gradients[1, k] <- -(0.5) * sum(diag(mat_temp))
# #   }
# #
# #   # Partial derivatives for period t >= 2
# #   for (i in 2:nrow(r)) {
# #     dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% r[(i-1), ] %*% t(r[(i-1), ])) + t(kronecker(g, g)) %*% dHda
# #     dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% ht) + t(kronecker(g, g)) %*% dHdg
# #     dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination) + t(kronecker(g, g)) %*% dHdc
# #
# #     dHdtheta <- t(cbind(dHdc, dHda, dHdg))
# #
# #     ht <- c_full + t(a) %*% r[(i-1), ] %*% t(r[(i-1), ]) %*% a + t(g) %*% ht %*%g
# #
# #     ht_sqrt_inv <- solve(sqrtm(ht))
# #     et <- ht_sqrt_inv %*% r[i,]
# #
# #     for (k in 1:nrow(theta)) {
# #       dh <- matrix(dHdtheta[k, ], N, N)
# #
# #       mat_temp <- ht_sqrt_inv %*% dh %*% ht_sqrt_inv %*% (diag(N) - tcrossprod(et))
# #
# #       gradients[i, k] <- -(0.5) * sum(diag(mat_temp))
# #     }
# #   }
# #
# #   return(gradients)
# # }
#
#



# Log-Likelihood function
# loglike_bekk <- function(theta, r) {
#   # convert to matrices
#   n <- ncol(r)
#   #Length of each series
#   NoOBs <- nrow(r)
#   numb_of_vars <- 2*n^2+n*(n+1)/2
#   C <- matrix(0,ncol = n,nrow = n)
#   index <- 1
#   for(i in 1:n){
#     for (j in i:n) {
#       C[j,i] <- theta[index]
#       index <- index+1
#     }
#   }
#   C <- t(C)
#   A = matrix(theta[index:(index+n^2-1)], n)
#   G = matrix(theta[(index+n^2):numb_of_vars], n)
#
#   # check constraints
#   if (!valid_bekk(C, A, G)) {
#     return(-1e25)
#   }
#
#   # compute inital H
#   H <- (t(r) %*% r) / nrow(r)
#
#   CC  <- t(C) %*% C
#   At  <- t(A)
#   Gt  <- t(G)
#
#
#   llv = numeric(NoOBs)
#   llv[1] <- -0.5 * (n * log(2 * pi) + log(det(H)) + c(t(r[1,]) %*% solve(H) %*% r[1,]))
#   for (i in 2:NoOBs) {
#     H <- (CC + At %*% (r[i - 1,] %*% t(r[i - 1,])) %*% A + Gt %*% H %*% G)
#     llv[i] = -0.5 * (n * log(2 * pi) + log(det(H)) + c(t(r[i,]) %*% solve(H) %*% r[i,]))
#   }
#   return(sum(llv))
# }

#Testing existance, uniqueness and stationarity
# valid_bekk <- function(C, A, G) {
#   # condition for positive-definit covariance
#   if (any(diag(C) < 0)) {
#     return(0)
#   }
#
#   # condition for uniqueness
#   if (A[1, 1] < 0 || G[1, 1] < 0) {
#     return(0)
#   }
#
#   # check stationarity for BEKK(1,1): Engle & Kroner (1995), Prop. 2.7
#   if (!all(abs(eigen(kronecker(A, A)
#                      + kronecker(G, G))$values) < 1)) {
#     return(0)
#   }
#
#   return(1)
# }

hessian_BEKK=function(theta,r){

  n=nrow(r)
  N <- ncol(r)
  N2 <- N^2
  NoOfVars_C = N*(N+1)/2

L_elimination <- elimination_mat(N)
D_duplication <- duplication_mat(N)
K_commutation <- commutation_mat(N)
D_gen_inv <- solve(crossprod(D_duplication)) %*% t(D_duplication)


C1=kronecker(kronecker(diag(N),K_commutation),diag(N))
C2=2*(kronecker(diag(N2),D_duplication%*%D_gen_inv))
C3=C2%*%C1

#BEKK-Matrices
c_mat <-  coef_mat(theta, 2)

c0 <- c_mat$c0
a <- c_mat$a
g <- c_mat$g

c_full=crossprod(c0)

#often done calculations
t_kron_g=t(kronecker(g,g))
gt=t(g)
kron_comm_cgt=kronecker(K_commutation,c(gt))


#Hessian
hessian = matrix(0, ncol = nrow(theta),nrow=nrow(theta))

#Second derivatives for t=1
dHdada = matrix(0,nrow=N2^2,ncol=N2) #Fehler bei h2??
dHdadc = matrix(0,nrow=N2^2,ncol=NoOfVars_C)
dHdadg=matrix(0,nrow=N2^2,ncol=N2)

dHdgda = matrix(0,nrow=N2^2,ncol=N2)
dHdgdg = matrix(0,nrow=N2^2,ncol=N2) #Fehler bei h2??
dHdgdc = matrix(0,nrow=N2^2,ncol=NoOfVars_C)

dHdcdc = matrix(0,nrow=NoOfVars_C*N2,ncol=NoOfVars_C)
dHdcdg = matrix(0,nrow=NoOfVars_C*N2,ncol=N2)
dHdcda = matrix(0,nrow=NoOfVars_C*N2,ncol=N2)

# Partial derivatives for initial period t = 1
ht <- crossprod(r)/nrow(r)

dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% ht)
dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(g) %*% ht)
dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination)

dHdtheta <- t(cbind(dHdc, dHda, dHdg))


ht_inv=solve(ht)

dHHdc=cbind(dHdcdc,dHdcda,dHdcdg)
dHHda=cbind(dHdadc,dHdada,dHdadg)
dHHdg=cbind(dHdgdc,dHdgda,dHdgdg)

dHHdtheta=rbind(dHHdc,dHHda,dHHdg)


matt=matrix(0, nrow=nrow(theta),ncol=N2)
matt[1,1]=1

for ( i in 2:nrow(theta)){
mat=matrix(0,nrow=nrow(theta),ncol=N2)
mat[i,1]=1
matt=cbind(matt,mat)
}
dHH=matt%*%dHHdtheta

for (i in 2:N2){
  matt=cbind(matrix(0,nrow=nrow(mat),ncol = 1),matt[,1:ncol(matt)-1])
  dHH=rbind(dHH,matt%*%dHHdtheta)
}

for (i in 1:nrow(theta)){
  for (j in 1:nrow(theta)){
    dhi = t(dHdtheta[i,])
    dhi = matrix(dhi,ncol=N,nrow=N);

    dhj = t(dHdtheta[j,])
    dhj = matrix(dhj,ncol=N,nrow=N);

    temp = dHH[i,j];

    for (k in 2:N2) {
    temp = cbind(temp, dHH[(k-1)*nrow(theta)+i,j])
    }

    temp = matrix(temp,ncol=N, nrow=N)
    ddh = t(temp)

    mat = ddh%*%ht_inv-dhi%*%ht_inv%*%dhj%*%ht_inv+r[1,]%*%t(r[1,])%*%ht_inv%*%dhj%*%ht_inv%*%dhi%*%ht_inv-
      r[1,]%*%t(r[1,])%*%ht_inv%*%ddh%*%ht_inv+r[1,]%*%t(r[1,])%*%ht_inv%*%dhi%*%ht_inv%*%dhj%*%ht_inv
    #update hessian
    hessian[i,j] = -0.5*sum(diag(mat))
  }
}

#Second derivatives for t>1

for (i in 2:n){
  dHdada = C3%*%(kronecker(c(diag(N)),kronecker(r[i-1,]%*%t(r[i-1,]),diag(N)))%*%K_commutation)+kronecker(diag(N2),t_kron_g) %*%dHdada
  dHdadg = kronecker(t(dHda), diag(N2))%*%C1%*%(kron_comm_cgt+kronecker(c(gt),K_commutation))+(kronecker(diag(N2),t(kronecker(g,gt))))%*%dHdadg
  dHdadc = dHdadc # Always zero (row may be deleted later on)

  dHdgda = C3%*%(kronecker(c(diag(N)),kronecker(diag(N),gt)%*%dHda ))+(kronecker(diag(N2),t_kron_g))%*%dHdgda
  dHdgdg = C3%*%(kronecker(c(diag(N)),kronecker(ht,diag(N))%*%K_commutation+kronecker(diag(N),gt)%*%dHdg ))+kronecker(t(dHdg),diag(N2))%*%C1%*%(kron_comm_cgt+kronecker(c(gt),K_commutation))+kronecker(diag(N2),t_kron_g)%*%dHdgdg
  dHdgdc = C3%*%kronecker(c(diag(N)),kronecker(diag(N),gt)%*%dHdc)+kronecker(diag(N2),t_kron_g)%*%dHdgdc

  dHdcda = dHdcda # Always zero (row may be deleted later on)
  dHdcdg = kronecker(t(dHdc),diag(N2))%*%C1%*%(kron_comm_cgt+kronecker(c(gt),K_commutation))+kronecker(diag(NoOfVars_C),t_kron_g)%*%dHdcdg
  dHdcdc = 2*kronecker(L_elimination,D_duplication%*%D_gen_inv)%*%C1%*%kronecker(diag(N2),c(diag(N)))%*%t(L_elimination)+kronecker(diag(NoOfVars_C),t_kron_g)%*%dHdcdc

  #updating first derivatives
  dHda <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), t(a) %*% r[(i-1), ] %*% t(r[(i-1), ])) + t(kronecker(g, g)) %*% dHda
  dHdg <- 2 * D_duplication %*% D_gen_inv %*% kronecker(diag(N), gt %*% ht) + t(kronecker(g, g)) %*% dHdg
  dHdc <- 2 * D_duplication %*% D_gen_inv %*% kronecker(t(c0), diag(N)) %*% t(L_elimination) + t(kronecker(g, g)) %*% dHdc

  #update h
  ht <- c_full + t(a) %*% r[(i-1), ] %*% t(r[(i-1), ]) %*% a + gt %*% ht %*%g

  #Coputation of Hessian for each t>1


  dHdtheta <- t(cbind(dHdc, dHda, dHdg))


  ht_inv=solve(ht)

  dHHdc=cbind(dHdcdc,dHdcda,dHdcdg)
  dHHda=cbind(dHdadc,dHdada,dHdadg)
  dHHdg=cbind(dHdgdc,dHdgda,dHdgdg)
  dHHdtheta=rbind(dHHdc,dHHda,dHHdg)


  matt=matrix(0, nrow=nrow(theta),ncol=N2)
  matt[1,1]=1

  for ( j in 2:nrow(theta)){
    mat=matrix(0,nrow=nrow(theta),ncol=N2)
    mat[j,1]=1
    matt=cbind(matt,mat)
  }
  dHH=matt%*%dHHdtheta

  for (j in 2:N2){
    matt=cbind(matrix(0,nrow=nrow(mat),ncol = 1),matt[,1:ncol(matt)-1])
    dHH=rbind(dHH,matt%*%dHHdtheta)
  }



  for (l in 1:nrow(theta)){
    for (j in 1:nrow(theta)){
      dhi = t(dHdtheta[l,])
      dhi = matrix(dhi,ncol=N,nrow=N);

      dhj = t(dHdtheta[j,])
      dhj = matrix(dhj,ncol=N,nrow=N);

      temp = dHH[l,j];

      for (k in 2:N2) {
        temp = cbind(temp, dHH[(k-1)*nrow(theta)+l,j])
      }

      temp = matrix(temp,ncol=N, nrow=N)
      ddh = t(temp)
      #get partial derivtatives for ll
      mat = ddh%*%ht_inv-dhi%*%ht_inv%*%dhj%*%ht_inv+r[i,]%*%t(r[i,])%*%ht_inv%*%dhj%*%ht_inv%*%dhi%*%ht_inv-
        r[i,]%*%t(r[i,])%*%ht_inv%*%ddh%*%ht_inv+r[i,]%*%t(r[i,])%*%ht_inv%*%dhi%*%ht_inv%*%dhj%*%ht_inv
      #update hessian
      hessian[l,j] = hessian[l,j]-0.5*sum(diag(mat))

    }
  }

}

return(hessian*(-1))
}

# Obtaining QML t-ratios
QML_t_ratios <- function(theta, r) {
  s1 <- score_bekk(theta, r)
  s1 <- crossprod(s1)

  s2 <- hessian_BEKK(theta, r)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(diag(s2))

  return(theta/s2)
}
