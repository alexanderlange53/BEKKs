test_that("Log-likelihood works with 2-dimensional sample data", {

  theta <-  c(0.026486941, -0.054371863, 0.118496142, 0.199628848, -0.007720701, -0.008865618,
              0.298887218, 0.976954018, 0.003170047, 0.019283272, 0.941393723)
  theta <- matrix(theta, nrow = 11)


  expect_equal(trunc(loglike_bekk(theta, data.matrix(StocksBonds))), -7425)
})

test_that("Log-likelihood works with 3-dimensional sample data", {

  theta <-  c( 0.0010122220, -0.0002678085,  0.0001905036,  0.0008937830,  0.0003111280,  0.0003893194,
               0.2453276126,  0.0222127631, -0.0434456548, -0.0366629406,  0.2637613982, -0.0248105673,
               0.0624407934, -0.0016446789,  0.1733957064,  0.9651528236, -0.0060973318,  0.0155928988,
               0.0125157899,  0.9606411455,  0.0008724460, -0.0196494968, -0.0009070321,  0.9801993180)
  theta <- matrix(theta, nrow = 24)


  expect_equal(trunc(loglike_bekk(theta, GoldStocksBonds)), 75249)
})


