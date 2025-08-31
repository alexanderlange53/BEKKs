test_that("Valid BEKK model 2-dims", {
  signs <- as.matrix(c(-1,-1))
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  B <- matrix(c(0.062342428, -0.000231618, 0.000453621,  0.0023456), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_true(valid_asymm_bekk(C,A,B,G, StocksBonds, signs))
})

test_that("Invalid BEKK model due to asymmetrics", {
  signs <- as.matrix(c(-1,-1))
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  B <- matrix(c(0.262342428, -0.000231618, 0.000453621,  0.0023456), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_asymm_bekk(C,A,B,G, StocksBonds, signs))
})

test_that("Invalid BEKK model due to non-uniquness 2-dims", {
  signs <- as.matrix(c(-1,-1))
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(-1.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  B <- matrix(c(0.062342428, -0.000231618, 0.000453621,  0.0023456), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_asymm_bekk(C,A,B,G, StocksBonds, signs))
})

test_that("Invalid BEKK model due to non stationarity 2-dims", {
  signs <- as.matrix(c(-1,1))
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -1.108865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  B <- matrix(c(0.462342428, -0.000231618, 0.000453621,  0.623456), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 1.51928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_asymm_bekk(C,A,B,G, StocksBonds, signs))
})
