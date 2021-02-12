
#H2 Code
gridSearch_BEKK=function(r){
  n=ncol(r)

  uncond_var=t(r)%*%r/nrow(r)
  a=matrix(0,ncol = n,nrow = n)
  g=matrix(0,ncol = n,nrow = n)
  c=matrix(0,ncol = n,nrow = n)
  #th0=numeric(2*n^2+n*(n+1)/2)

  for (i in 1:n){
    a[i,i] = sqrt(0.05)
    g[i,i] = sqrt(0.9)
    c[i,i] = 0.05*uncond_var[i,i]
  }

  for (i in 1:n){
    for (j in seq(i,n)){

      cij = uncond_var[i,j]/sqrt(uncond_var[i,i]*uncond_var[j,j])
      c[i,j] = cij*sqrt(c[i,i]*c[j,j])
      c[j,i] = c[i,j]

    }
  }

  c = t(chol(c))
  c0 = c[,1]

for (i in 2:(n-1)){
  c0 = c(c0,c[i:n,i])
}

c0 = c(c0,c[n,n])

# deterministic variance components

th0 = c(c0,c(a),c(g))

#change elements of a and g and compute likelihood in each step

likmax = -9999999

#print th0
result= recursive_search_BEKK(r,c0,c(a),c(g),1,th0,likmax)
th0=result[[1]]
likmax=result[[2]]
return(list(th0,likmax))

}


recursive_search_BEKK=function(r,c0,avec,gvec,index,thetaopt,likmax){

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
    avec[index,1] = val
  } else{
    gvec[index-n^2,1] = val
  }
#last element is excluded
  if (index < (2*n^2-1)){
# recursive step
 result= recursive_search_BEKK(r,c0,avec,gvec,index+1,thetaopt,likmax)
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
#     for (j in 1:numb_of_vars){
#       if(j==counter & j<=(n*(n+1)/2)){
#         theta[j]=runif(1,min = 0, max = 1)
#         counter=counter+diagonal_elements
#         diagonal_elements=diagonal_elements-1
#
#       }else if(j==(n*(n+1)/2+1+diagonal_counter*n)){
#         theta[j]=runif(1,min = 0, max = 1)
#
#         diagonal_counter=diagonal_counter+1
#       }else{
#         theta[j]=runif(1,min = -0.5, max = 0.5)
#       }
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
#   return(list(thetaOptim,best_val))
# }
