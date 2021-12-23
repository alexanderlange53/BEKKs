test_that("Testing GridSearch works with StocksBonds", {
  set.seed(6723217)
  gs1 <- random_grid_search_BEKK(data.matrix(StocksBonds))

  # Starting values should be equal for same seed

  expect_equal(trunc(sum(gs1[[1]])), 0 )
})
