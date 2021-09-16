#' Estimating a BEKK(1, 1) model
#'
#' @param r data input
#' @param init_values initial values for BEKK parameter
#' @param QML_t_ratios Logical. If QML_t_ratios = 'TRUE', the t-ratios of the BEKK parameter matrices
#'                     are exactly calculated via second order derivatives.
#' @param max_iter maximum number of BHHH algorithm iterations
#' @param crit determiens the precision of the BHHH algorithm
#'
#' @examples
#' \donttest{
#'
#' data(bivariate)
#' x1 <- bekk(BI, init_values = NULL,
#' QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#'
#' summary(x1)
#'
#' plot(x1)
#'
#' }
#' @export

bekk <- function(r, init_values = NULL, QML_t_ratios = FALSE,
                 seed = NULL, max_iter = 50, crit = 1e-9, nc = 1, asymmetric=FALSE){

  # Checking for valid input
  r <- as.matrix(r)
  if (any(is.na(r))) {
    stop("\nNAs in y.\n")
  }
  if (ncol(r) < 2) {
    stop("The matrix 'r' should contain at least two variables.")
  }
  if (is.null(colnames(r))) {
    colnames(r) <- paste("y", 1:ncol(r), sep = "")
  }


  N <- ncol(r)
if(asymmetric==FALSE){
  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_BEKK(r)
      theta <- theta[[1]]
  } else if (init_values == 'random') {
      if(is.null(seed)) {
        seed <- round(runif(nc, 1, 100))
      } else {
        set.seed(seed)
        seed <- round(runif(nc, 1, 100))
      }
    cat('Generating starting values \n')
      theta_list <- pblapply(seed, random_grid_search_BEKK, r = r,
                             nc = nc,
                             cl =nc)
      max_index <- which.max(sapply(theta_list, '[[', 'best_val'))
      theta <- theta_list[[max_index]]
      theta <- theta[[1]]
  } else if (init_values == 'simple') {
    uncond_var <- crossprod(r)/nrow(r)
    A <- matrix(0, ncol = N, nrow = N)
    G <- matrix(0, ncol = N, nrow = N)
    C <- matrix(0, ncol = N, nrow = N)
    #th0=numeric(2*n^2+n*(n+1)/2)

    diag(A) <- 0.3
    diag(G) <- 0.92
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

    if (N > 2) {
      for (i in 2:(N-1)){
        C0 = c(C0, C[i:N, i])
      }
    }

    C0 = c(C0, C[N, N])

    theta = c(C0, c(A), c(G))

    }
  } else {
    if(length(init_values) != 2 * N^2 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_bekk(r, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios(params$theta, r)
    tratios_mat <- coef_mat(tratios, N)
  } else {
    tratios_mat <- coef_mat(params$t_val, N)
  }

  param_mat <- coef_mat(params$theta, N)

  var_process <- sigma_bekk(r, param_mat$c0, param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional standard deviation of \n', colnames(r)[k])
        #sigma_t[,k2] <- sqrt(sigma_t[,k2])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of \n', colnames(r)[i], ' and ', colnames(r)[j])
        #sigma_t[,k2] <- sigma_t[,k2]/(sigma_t[,i*(k-1)] * sqrt(sigma_t[,j*N]))
        k2 <- k2 +1
      }
    }
  }

  for (i in 1:nrow(sigma_t)) {
    tm <- matrix(unlist(sigma_t[i,]), N, N, byrow = T)
    tm2 <- sqrt(solve(diag(diag(tm))))%*%tm%*%sqrt(solve(diag(diag(tm))))
    diag(tm2) <- sqrt(diag(tm))
    sigma_t[i,] <- c(tm2)
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]


  if (inherits(r, "ts")) {
    sigma_t <- ts(sigma_t, start = time(r)[1], frequency = frequency(r))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_bekk(param_mat$c0, param_mat$a, param_mat$g)


  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 G_t = tratios_mat$g,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 e_t = var_process$e_t,
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 data = r)
  class(result) <- 'bekk'
  return(result)
}else if(asymmetric==TRUE) {
  if(!is.numeric(init_values)) {
    if (is.null(init_values)) {
      theta <- gridSearch_asymmetricBEKK(r)
      theta <- theta[[1]]
    } else if (init_values == 'random' && nc>1) {
      if(is.null(seed) ) {
        seed <- round(runif(nc, 1, 100))
      } else {
        set.seed(seed)
        seed <- round(runif(nc, 1, 100))
      }

      cat('Generating starting values \n')
      theta_list <- pblapply(X=seed, FUN=random_grid_search_BEKK, r = r,
                             nc = nc,cl=cl         )
      max_index <- which.max(sapply(theta_list, '[[', 'best_val'))
      theta <- theta_list[[max_index]]
      theta <- theta[[1]]
    } else if (init_values == 'random' && nc==1) {
      if(is.null(seed) ) {
        seed <- round(runif(1, 1, 100))
      } else {
        set.seed(seed)
        seed <- round(runif(1, 1, 100))
      }

      cat('Generating starting values \n')
      theta_max <- random_grid_search_BEKK(r=r,seed=seed,nc = nc)

      theta=theta_max$thetaOptim

    } else if (init_values == 'simple') {
      uncond_var <- crossprod(r)/nrow(r)
      A <- matrix(0, ncol = N, nrow = N)
      B <- matrix(0, ncol = N, nrow = N)
      G <- matrix(0, ncol = N, nrow = N)
      C <- matrix(0, ncol = N, nrow = N)
      #th0=numeric(2*n^2+n*(n+1)/2)

      diag(A) <- 0.25
      diag(B) <- 0.05
      diag(G) <- 0.92
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

      if (N > 2) {
        for (i in 2:(N-1)){
          C0 = c(C0, C[i:N, i])
        }
      }

      C0 = c(C0, C[N, N])

      theta = c(C0, c(A), c(G))

    }
  } else {
    if(length(init_values) != 3 * N^2 + N * (N + 1)/2) {
      stop('Number of initial parameter does not match dimension of data.')
    }
    theta <- init_values
  }

  theta <- matrix(theta, ncol =1)

  params <- bhh_asymm_bekk(r, theta, max_iter, crit)

  if (QML_t_ratios == TRUE) {
    tratios <- QML_t_ratios(params$theta, r)
    tratios_mat <- coef_mat_asymm(tratios, N)
  } else {
    tratios_mat <- coef_mat_asymm(params$t_val, N)
  }

  param_mat <- coef_mat_asymm(params$theta, N)

  var_process <- sigma_bekk(r, param_mat$c0, param_mat$a, param_mat$g)
  sigma_t <- as.data.frame(var_process$sigma_t)
  colnames(sigma_t) <- rep(1, N^2)

  k <- 1
  k2 <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        colnames(sigma_t)[k2] <- paste('Conditional SD of', colnames(r)[k])
        k <- k + 1
        k2 <- k2 +1
      } else {
        colnames(sigma_t)[k2] <- paste('Conditional correlation of', colnames(r)[i], ' and ', colnames(r)[j])
        k2 <- k2 +1
      }
    }
  }

  elim <- elimination_mat(N)
  sigma_t <- sigma_t[, which(colSums(elim) == 1)]

  if (inherits(r, "ts")) {
    sigma_t <- ts(sigma_t, start = time(r)[1], frequency = frequency(r))
  }

  # Final check if BEKK is valid
  BEKK_valid <- valid_asymm_bekk(param_mat$c0, param_mat$a, param_mat$b, param_mat$g)


  params$likelihood_iter <- params$likelihood_iter[params$likelihood_iter != 0]

  result <- list(C0 =  param_mat$c0,
                 A = param_mat$a,
                 B = param_mat$b,
                 G = param_mat$g,
                 C0_t = tratios_mat$c0,
                 A_t = tratios_mat$a,
                 B_t = tratios_mat$b,
                 G_t = tratios_mat$g,
                 log_likelihood = params$likelihood,
                 BEKK_valid = BEKK_valid,
                 sigma_t = sigma_t,
                 e_t = var_process$e_t,
                 iter = params$iter,
                 likelihood_iter = params$likelihood_iter,
                 data = r)
  class(result) <- 'asymmetricbekk'
  return(result)
}
}
