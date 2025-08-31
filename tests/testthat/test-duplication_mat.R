test_that("duplication matrix n = 2", {
  D <- duplication_mat(2)

  m <- matrix(c(1, 0, 0,
                0, 1, 0,
                0, 1, 0,
                0, 0, 1), ncol = 3, byrow = TRUE)

  expect_equal(D, m)
})

test_that("duplication matrix n = 3", {
  D <- duplication_mat(3)

  m <- matrix(c(1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1), ncol = 6, byrow = TRUE)

  expect_equal(D, m)
})
