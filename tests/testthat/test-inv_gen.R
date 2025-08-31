test_that("inv_gen works with nonsingular positive definite 2x2 matrix", {
  testm <- matrix(c(0.07188, -0.0124, -0.0124, 0.833),2,2)

  expect_equal(inv_gen(testm), solve(testm))
})

test_that("inv_gen works with singular not positive definite 2x2 matrix", {
  testm <- matrix(c(0.07188, -0.01241, 0, 0),2,2)

  expect_equal(round(inv_gen(testm), 2), matrix(c(13.51, 0, -2.331, 0), 2, 2), tolerance =0.1)
})

