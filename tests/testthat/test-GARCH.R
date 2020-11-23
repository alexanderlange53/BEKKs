test_that("Testing GARCH(1,1) works", {
  # Replicating results from book statistics of financial markets 2019 page 290, table 13.4

  g1 <- garch(DAX30, 1, 1)

  # Likelihood (book says 7614.9)
  expect_equal(round(g1$loglik, 1), 7617.4)
  # Parameter omega, alpha and beta
  expect_equal(round(g1$theta, 3),  c(0, 0.079, 0.914))
  # t-values
  expect_equal(round(g1$t_values, 1),  c(4.3, 10.0, 112.6))
})
