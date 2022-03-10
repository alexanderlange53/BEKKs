test_that("score_dbekk works with 2-dimensional test data set", {
  init <- c(0.062, -0.02689, 0.210177, 0.2236, 0.2236, 0.94868, 0.9487)

  theta <- matrix(init, nrow = 7)

  score <- score_dbekk(theta, data.matrix(StocksBonds))
  score <- sum(score)
  expect_equal(trunc(score), -21797)
})

test_that("score_dbekk works with 3-dimensional test data set", {
  theta=c(0.0015876334,  0.0002258635,  0.0006750149,  0.0014719950, -0.0005589923,  0.0012654071,  0.2610727800, 0.2989689227,
          0.2682630137,  0.9561381437,  0.9440457439, 0.9599385388)
  score <- score_dbekk(as.matrix(theta), GoldStocksBonds)
  score <- sum(score)
  expect_equal(trunc(score), -1054355)
})
