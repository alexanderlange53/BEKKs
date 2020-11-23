test_that("Testing GARCH(1,1) works with DAX data", {
  # Replicating results from book statistics of financial markets 2019 page 290, table 13.4

  g1 <- garch(DAX30, 1, 1)

  # Likelihood (book says 7614.9)
  expect_equal(round(g1$loglik, 1), 7617.4)
  # Parameter omega, alpha and beta
  expect_equal(round(g1$theta, 3),  c(0, 0.079, 0.914))
  # t-values
  expect_equal(round(g1$t_values, 1),  c(4.3, 10.0, 112.6))
})

test_that("Testing GARCH(2,1) works with DAX data", {

  g2 <- garch(DAX30, 2, 1)

  # Likelihood (book says 7614.9)
  expect_equal(round(g2$loglik, 1), 7584)
  # Parameter omega, alpha and beta
  expect_equal(round(g2$theta, 3),  c(0, 0.079, 0.914))
  # t-values
  expect_equal(round(g2$t_values, 1),  c(4.3, 10.0, 112.6))
})
