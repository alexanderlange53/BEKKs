# Funtion to extract conditional standard deviations from bekk_fit object
extract_csd <- function(x) {

  csd <- matrix(NA, nrow = nrow(x$sigma_t), ncol = ncol(x$data))
  csd_names <- rep(NA, ncol(x$data))

  counter <- 1
  for (i in 1:ncol(x$data)) {
    csd[, i] <- x$sigma_t[, counter]
    csd_names[i] <- colnames(x$sigma_t)[counter]
    counter <- counter + ncol(x$data) - (i - 1)
  }

  colnames(csd) <- csd_names
  return(csd)
}
