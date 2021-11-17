bekk_mc_eval <- function(object, spec, sample_sizes, iter, nc = 1) {
  UseMethod('bekk_mc_eval')
}

bekk_mc_eval.bekk <- function(object, spec, sample_sizes, iter, nc = 1) {
    xx <- process_object(object)
    theta <- xx$theta

    mse <- rep(NA, length(sample_sizes))
    index <- 1

    for(j in sample_sizes) {
      print(paste('Sample size: ', j))
      sim_dat <- vector(mode = "list", iter)

      sim_dat <- mclapply(1:iter, function(x){bekk_sim(object, nobs = j)},  mc.cores = nc)

      dd <- pblapply(sim_dat, function(x){bekk_fit(spec = spec, data = x)}, cl = nc)

      mse[index] <- sum(unlist(lapply(dd, RMSE, theta_true = theta)))/iter
      index <- index +1

    }

    result <- data.frame(Sample_size = sample_sizes, MSE = mse)
    colnames(result) <- c('Sample', 'MSE')
    result <- list(result)

    class(result) <- 'bekkMC'
    return(result)
}

RMSE <- function(x, theta_true) {
  theta_est <- x$theta
  return(sum(sqrt(((theta_true - theta_est) / theta_true)^2)))
}

plot.bekkMC <- function(x, ...) {
  msep <- x[[1]]
  ggplot(msep) + geom_line(aes(x = Sample, y = MSE)) + theme_bw()
}
