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

  expect_equal(round(g2$loglik, 1), 7615.3)
  # Parameter omega, alpha and beta
  expect_equal(round(g2$theta, 3),  c(0.000, 0.063, 0.922, 0.007))
  # t-values
  expect_equal(round(g2$t_values, 1),  c(2.8, 3.1, 2.7, 0.0))
})

test_that("Testing GARCH(1,2) works with DAX data", {

  g3 <- garch(DAX30, 1, 2)

  expect_equal(round(g3$loglik, 1), 7620.7)
  # Parameter omega, alpha and beta
  expect_equal(round(g3$theta, 3),  c(0.000, 0.014, 0.080, 0.895))
  # t-values
  expect_equal(round(g3$t_values, 1),  c(4.7, 1.1, 4.8, 86.6))
})

test_that("Testing GARCH(2,2) works with DAX data", {

  g4 <- garch(DAX30, 2, 2)

  expect_equal(round(g4$loglik, 1), 7617.2)
  # Parameter omega, alpha and beta
  expect_equal(round(g4$theta, 3),  c(0.000, 0.012, 0.059, 0.917, 0.002))
  # t-values
  expect_equal(round(g4$t_values, 1),  c(2.3, 0.8, 1.9, 2.1, 0.0))
})

test_that("Testing GARCH(3,1) works with DAX data", {

  g5 <- garch(DAX30, 3, 1)

  expect_equal(round(g5$loglik, 1), 7584.4)
  # Parameter omega, alpha and beta
  expect_equal(round(g5$theta, 3),  c(0.000, 0.048, 0.701, 0.229, 0.001))
  # t-values
  expect_equal(round(g5$t_values, 1),  c(2.3, 2.4, 1.6, 0.4, 0.0))
})

test_that("Testing GARCH(1,3) works with DAX data", {

  g6 <- garch(DAX30, 1, 3)

  expect_equal(round(g6$loglik, 1), 7618.5)
  # Parameter omega, alpha and beta
  expect_equal(round(g6$theta, 3),  c(0.000, 0.012, 0.067, 0.021, 0.890))
  # t-values
  expect_equal(round(g6$t_values, 1),  c(4.6, 0.9, 2.5, 0.8, 73.4))
})

test_that("Testing GARCH(3,3) works with DAX data", {

  g7 <- garch(DAX30, 3, 3)

  expect_equal(round(g7$loglik, 1), 7521.6)
  # Parameter omega, alpha and beta
  expect_equal(round(g7$theta, 3),  c(0.000, 0.014, 0.033, 0.023, 0.005, 0.562, 0.310))
  # t-values
  expect_equal(round(g7$t_values, 1),  c(1.0, 0.7, 1.0, 0.5, 0.0, 0.4, 0.3))
})
