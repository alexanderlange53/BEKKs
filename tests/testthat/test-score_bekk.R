test_that("score_bekk works with 2-dimensional test data set", {
  init <- c(0.062, -0.02689, 0.210177, 0.2236, -0.03, -0.03,
            0.2236, 0.94868, 0.03, -0.03, 0.9487)

  theta <- matrix(init, nrow = 11)

  score <- score_bekk(theta, data.matrix(StocksBonds))
  score <- sum(score)
  expect_equal(trunc(score), -166954)
})
