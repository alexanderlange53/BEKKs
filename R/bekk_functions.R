process_object <- function(x) {
  UseMethod('process_object')
}

process_object.bekkFit <- function(x) {
  theta <- x$theta
  N <- ncol(x$C0)
  signs <- x$signs
  BEKK_valid <- x$BEKK_valid
  expected_signs <- x$expected_signs

  return(list(theta = theta,
              N = N, signs=signs, expected_signs = expected_signs, BEKK_valid = BEKK_valid
              ))
}

process_object.bekkSpec <- function(x) {
  if (is.null(x$init_values)) {
    stop('Please provide "initial_values" in "bekk_spec" as paramater for simulation.')
  }
  if (is.null(x$N)) {
    stop('Please provide "N" in "bekk_spec" as dimension for simulation.')
  }

  theta <- x$init_values
  N <- x$N

  if(is.null(x$signs) && x$model$asymmetric == TRUE){
    signs=as.matrix(rep(-1,N))

  }else{
    signs=x$signs
  }

  if(x$model$asymmetric == FALSE){
    par = coef_mat(theta,N)
    BEKK_valid = valid_bekk(par$c0, par$a, par$g)
  } else{
    par = coef_mat_asymm(theta,N)
    BEKK_valid = valid_asymm_bekk_sim(par$c0, par$a, par$b, par$g, 1/(N^2),signs)
  }

  expected_signs=1/(N^2)
  return(list(theta = theta,
              N = N, signs = signs, expected_signs = expected_signs, BEKK_valid = BEKK_valid))
}

# Obtaining QML t-ratios
QML_t_ratios <- function(theta, r) {
  s1 <- score_bekk(theta, r)
  s1 <- crossprod(s1)

  s2 <- hesse_bekk(theta, r)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

QML_t_ratios_asymm <- function(theta, r, signs) {
  s1 <- score_asymm_bekk(theta, r, signs)
  s1 <- crossprod(s1)

  s2 <- hesse_asymm_bekk(theta, r, signs)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

QML_t_ratios_dbekk <- function(theta, r) {
  s1 <- score_dbekk(theta, r)
  s1 <- crossprod(s1)

  s2 <- hesse_dbekk(theta, r)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

QML_t_ratios_dbekk_asymm <- function(theta, r, signs) {
  s1 <- score_asymm_dbekk(theta, r, signs)
  s1 <- crossprod(s1)

  s2 <- hesse_asymm_dbekk(theta, r, signs)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

QML_t_ratios_sbekk <- function(theta, r) {
  s1 <- score_sbekk(theta, r)
  s1 <- crossprod(s1)

  s2 <- hesse_sbekk(theta, r)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

QML_t_ratios_sbekk_asymm <- function(theta, r, signs) {
  s1 <- score_asymm_sbekk(theta, r, signs)
  s1 <- crossprod(s1)

  s2 <- hesse_asymm_sbekk(theta, r, signs)
  s2 <- solve(s2) %*% s1 %*% solve(s2)

  s2 <- sqrt(abs(diag(s2)))

  return(abs(theta/s2))
}

