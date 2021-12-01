test_that("Log-likelihood works with 2-dimensional sample data", {

  theta <-  c(0.026486941, -0.054371863, 0.118496142, 0.199628848, -0.007720701, -0.008865618,
              0.298887218, 0.976954018, 0.003170047, 0.019283272, 0.941393723)
  theta <- matrix(theta, nrow = 11)


  expect_equal(trunc(loglike_bekk(theta, data.matrix(StocksBonds))), -7425)
})

