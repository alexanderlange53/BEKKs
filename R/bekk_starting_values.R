
#H2 Code
gridSearch_BEKK <- function(r){
  N  <- ncol(r)

  uncond_var <- crossprod(r)/nrow(r)
  A <- matrix(0, ncol = N, nrow = N)
  G <- matrix(0, ncol = N, nrow = N)
  C <- matrix(0, ncol = N, nrow = N)
  #th0=numeric(2*n^2+n*(n+1)/2)

  diag(A) <- sqrt(0.05)
  diag(G) <- sqrt(0.9)
  diag(C) <- 0.05*diag(uncond_var)


  for (i in 1:N){
    for (j in seq(i,N)){

      cij <- uncond_var[i, j]/sqrt(uncond_var[i, i]*uncond_var[j, j])
      C[i,j] <- cij*sqrt(C[i, i]*C[j, j])
      C[j,i] <- C[i, j]

    }
  }

  C = t(chol(C))
  C0 = C[,1]

  if ((N-1) > 2) {
    for (i in 2:(N-1)){
      C0 = c(c0, C[i:N, i])
    }
  }

  C0 = c(C0, C[N, N])

  # deterministic variance components

  th0 = c(C0, c(A), c(G))

  # change elements of A and G and compute likelihood in each step

  likmax = -1e25

  #print th0
  result= recursive_search_BEKK(r, C0, c(A), c(G), 1, th0, likmax)
  th0=result[[1]]
  likmax=result[[2]]
  return(list(th0,likmax))

}


recursive_search_BEKK=function(r, c0, avec, gvec, index, thetaopt, likmax){

n=ncol(r)
start= -3
endr = 3
step = 6
indextest=0

if (index == n^2){
  index = index + 2;
} else if (index < n^2){
  indextest = (index-1)/(n+1)
} else{
  indextest = (index-n^2-1)/(n+1)
}
# we have a diagonal element
if (indextest - floor(indextest) == 0){
    index = index + 1
}

for (i in seq(start, endr, step)){
  val = i/100

  #set a and g respectively according to index, exclude diagonal elements
  if (index <= n^2){
    avec[index] = val
  } else{
    gvec[index-n^2] = val
  }
  #last element is excluded
  if (index < (2*n^2-1)){
    # recursive step
    result= recursive_search_BEKK(r, c0, avec, gvec, index+1, thetaopt, likmax)
    thetaopt=result[[1]]
    likmax=result[[2]]
  } else{
    #final step
    theta = c(c0,avec,gvec)
    #likelihood
    lik = loglike_bekk(theta,r)
    if (lik > likmax){
      thetaopt = theta
      likmax = lik
    }

  }
}

return(list(thetaopt,likmax))

}



#Now Grid Search with generation of random candidates
# random_grid_search_BEKK=function(r,sampleSize){
#   n=ncol(r)
#   numb_of_vars=2*n^2+n*(n+1)/2
#   theta=numeric(numb_of_vars)
#   thetaOptim=theta
#   best_val=-999999
#   #Generating random values for A, C and G
#
#   for (i in 1:sampleSize){
#     counter=1
#     diagonal_elements=n
#     diagonal_counter=0
#     for (j in 1:(n*(n+1)/2)){
#       if(j==counter & j<=(n*(n+1)/2)){
#         theta[j]=runif(1,min = 0, max = 1)
#         counter=counter+diagonal_elements
#         diagonal_elements=diagonal_elements-1
#
#       }else{
#         theta[j]=runif(1,min = -0.9, max = 0.9)
#       }
#     }
#     for(j in (1+n*(n+1)/2):numb_of_vars){
#       if(j==1+n*(n+1)/2 || j==(1+n*(n+1)/2 +n*n)){
#         theta[j]=runif(1,min = 0, max = 0.9)
#       }else{
#         theta[j]=runif(1,min = -0.9, max = 0.9)
#         }
#
#     }
#     C=matrix(0,ncol = n,nrow = n)
#     index=1
#     for(j in 1 : n){
#       for (k in j:n) {
#         C[k,j]=theta[index]
#         index=index+1
#       }
#     }
#
#
#     A = matrix(theta[index:(index+n^2-1)], n)
#     G = matrix(theta[(index+n^2):numb_of_vars], n)
#
#     if(!valid_bekk(C,A,G)){
#       NULL
#
#     }else{
#       llv=loglike_bekk(theta,r)
#       if(llv>best_val){
#         best_val=llv
#         thetaOptim=theta
#       }
#     }
#
#   }
#   return(list(thetaOptim,best_val,C,A,G))
# }

