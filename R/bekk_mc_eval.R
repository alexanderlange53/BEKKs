bekk_mc_eval <- function(spec, sample_sizes, iter, nc = 1) {
  if (class(spec) == 'bekk') {
    theta <- spec$theta
    N <- ncol(spec$C0)
    mse <- rep(NA, length(sample_sizes))
    index <- 1

    for(j in sample_sizes) {
      sim_dat <- vector(mode = "list", iter)

      # for(i in 1:iter) {
      #   sim_dat[[i]] <- bekk(bekk_sim(spec, j), nc =nc)
      # }

      for(i in 1:iter) {
        sim_dat[[i]] <- bekk_sim(spec, j)
      }

      dd <- pblapply(sim_dat, function(x){bekk(x)}, cl = nc)


      mse[index] <- sum(unlist(lapply(dd, RMSE, theta_true = theta)))/iter
      index <- index +1

    }

    result <- data.frame(Sample_size = sample_sizes, MSE = mse)
    colnames(result) <- c('Sample', 'MSE')
    result <- list(result)
    #rownames(result) <- sample_sizes

    class(result) <- 'bekkMC'
    return(result)

  } else if (class(spec) == 'bekkSimSpec') {

  } else {
    stop("The object 'spec' should be of class 'bekk' or 'bekkSimSpec'.")
  }
}

RMSE <- function(x, theta_true) {
  theta_est <- x$theta
  return(sum(sqrt(((theta_true - theta_est) / theta_true)^2)))
}

plot.bekkMC <- function(x, ...) {
  msep <- x[[1]]
  ggplot(msep) + geom_line(aes(x = Sample, y = MSE)) + theme_bw()
}
