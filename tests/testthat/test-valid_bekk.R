test_that("Valid BEKK model 2-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_true(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to negtaive heteroskedasticity 2-dims", {
  C <-  matrix(c(-0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to non uniquness 2-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(-1.199628848, -0.008865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to non stationarity 2-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0, 0.11849614), 2,2, byrow = T)
  A <- matrix(c(0.199628848, -1.108865618, -0.007720701,  0.298887218), 2, 2, byrow = T)
  G <- matrix(c(0.976954018, 1.51928327, 0.003170047, 0.94139372), 2, 2, byrow = T)

  expect_false(valid_bekk(C,A,G))
})

test_that("Valid BEKK model 3-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0.1, 0, 0.1,  0.11849614, 0, 0, 0.3), 3,3, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, 0.13, -0.007720701,  0.098887218, 0.1, -0.1, 0.1, 0.1), 3, 3, byrow = T)
  G <- matrix(c(0.76954018, 0.01928327, 0.003170047, 0.94139372, 0.4, -0.1, -0.01, 0.01, 0.2), 3, 3, byrow = T)

  expect_true(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to negtaive heteroskedasticity 2-dims", {
  C <-  matrix(c(-0.02648694, -0.05437186, 0.1, 0, 0.1,  0.11849614, 0, 0, 0.3), 3,3, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, 0.13, -0.007720701,  0.098887218, 0.1, -0.1, 0.1, 0.1), 3, 3, byrow = T)
  G <- matrix(c(0.76954018, 0.01928327, 0.003170047, 0.94139372, 0.4, -0.1, -0.01, 0.01, 0.2), 3, 3, byrow = T)

  expect_false(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to non uniquness 2-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0.1, 0, 0.1,  0.11849614, 0, 0, 0.3), 3,3, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, 0.13, -0.007720701,  0.098887218, 0.1, -0.1, 0.1, 0.1), 3, 3, byrow = T)
  G <- matrix(c(-1.76954018, 0.01928327, 0.003170047, 0.94139372, 0.4, -0.1, -0.01, 0.01, 0.2), 3, 3, byrow = T)

  expect_false(valid_bekk(C,A,G))
})

test_that("Invalid BEKK model due to non stationarity 2-dims", {
  C <-  matrix(c(0.02648694, -0.05437186, 0.4, 0, 0.1,  0.11849614, 0, 0, 0.3), 3,3, byrow = T)
  A <- matrix(c(0.199628848, -0.008865618, 0.33, -0.007720701,  0.298887218, 0.1, -0.3, 0.1, 0.11), 3, 3, byrow = T)
  G <- matrix(c(0.976954018, 0.01928327, 0.003170047, 0.94139372, 0.4, -0.1, -0.33, 0.01, 0.4), 3, 3, byrow = T)

  expect_false(valid_bekk(C,A,G))
})
