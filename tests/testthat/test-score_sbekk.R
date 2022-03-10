test_that("score_sbekk works with 2-dimensional test data set", {
  init <- c(0.062, -0.02689, 0.210177, 0.1346, 0.74868)
  theta <- matrix(init, nrow = 5)

  score <- score_sbekk(theta, data.matrix(StocksBonds))
  score <- sum(score)
  expect_equal(trunc(score), 55899)
})

test_that("score_sbekk works with 3-dimensional test data set", {

  init=c(0.0015876334,  0.0002258635,  0.0006750149,  0.0014719950, -0.0005589923,  0.0012654071,   0.082,  0.8599385388)
  theta <- matrix(init, nrow = 8)
  score <- score_sbekk(theta, GoldStocksBonds)
  score <- sum(score)
  expect_equal(trunc(score), 2000803)
})
