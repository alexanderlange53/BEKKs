test_that("bhh works with 2-dimensional test set", {
  theta <- c(0.062, -0.02689, 0.210177, 0.2236, -0.03, -0.03,
             0.2236, 0.94868, 0.03, -0.03, 0.9487)
  theta <- matrix(theta, nrow = 11)
  res <- bhh_bekk(data.matrix(StocksBonds), theta, max_iter = 50, crit = 1e-9)

  expect_equal(trunc(res$likelihood), -7382)
  expect_equal(res$iter, 13)
})
