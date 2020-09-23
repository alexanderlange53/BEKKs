test_that("ARCH(1) model with daily DAX returns works", {
  a1 <- arch(DAX30, 1)

  expect_equal(round(a1$loglik, 1), 7241.4)
  expect_equal(round(a1$theta, 3),  c(0, 0.295))
})

test_that("ARCH(3) model with daily DAX returns works", {
  a1 <- arch(DAX30, 3)

  expect_equal(round(a1$loglik, 1), 7433)
  expect_equal(round(a1$theta, 3),  c(0, 0.116, 0.267, 0.265))
})

test_that("ARCH(6) model with daily DAX returns works", {
  a1 <- arch(DAX30, 6)

  expect_equal(round(a1$loglik, 1), 7537)
  expect_equal(round(a1$theta, 3),  c(0.000, 0.038, 0.144, 0.163, 0.132, 0.175, 0.165))
})
