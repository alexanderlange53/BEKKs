test_that("score_asymm_bekk works with 2-dimensional test data set", {
  init <- c(0.062, -0.02689, 0.210177, 0.2236, -0.03, -0.03,
            0.2236, 0.05, -0.0005536, 0.00672,0.2753, 0.94868, 0.03, -0.03, 0.9487)
  signs = as.matrix(c(-1, -1))
  theta <- matrix(init, nrow = 15)

  score <- score_asymm_bekk(theta, data.matrix(StocksBonds),signs)
  score <- sum(score)
  expect_equal(trunc(score), -181125)
})

test_that("score_asymm_bekk works with 3-dimensional test data set", {
  signs = as.matrix(c(-1, -1, -1))

  init=c(0.0015876334,  0.0002258635,  0.0006750149,  0.0014719950, -0.0005589923,  0.0012654071,  0.2610727800,  0.0003034200,  0.0224087426,  0.0156283014,  0.2989689227,  0.0275966405,
          0.0131739627,  0.0495781047,  0.2682630137,
         0.12003, -0.0012, 0.0034234, 0.0066, 0.09234,  0.000234, -0.0000567, 0.00562146, 0.063466,
         0.9561381437,  0.0005135573, -0.0072066784,
          -0.0072261195,  0.9440457439, -0.0043198682, -0.0139102556, -0.0069184635,  0.9599385388)
  theta <- matrix(init, nrow = 33)
  score <- score_asymm_bekk(theta, GoldStocksBonds, signs)
  score <- sum(score)
  expect_equal(trunc(score), -805405)
})
