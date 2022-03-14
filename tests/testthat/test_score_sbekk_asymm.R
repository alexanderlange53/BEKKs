test_that("score_asymm_sbekk works with 2-dimensional test data set", {
  init <- c(0.062, -0.02689, 0.210177, 0.2236, 0.12, 0.94868)
  signs = as.matrix(c(-1, -1))
  theta <- matrix(init, nrow = 6)

  score <- score_asymm_sbekk(theta, data.matrix(StocksBonds),signs)
  score <- sum(score)
  expect_equal(trunc(score), -146558)
})

test_that("score_asymm_sbekk works with 3-dimensional test data set", {
  signs = as.matrix(c(-1, -1, -1))

  init=c(0.0015876334,  0.0002258635,  0.0006750149,  0.0014719950, -0.0005589923,  0.0012654071,  0.2610727800,  0.082,  0.9599385388)
  theta <- matrix(init, nrow = 9)
  score <- score_asymm_sbekk(theta, GoldStocksBonds, signs)
  score <- sum(score)
  expect_equal(trunc(score), -2393928)
})
