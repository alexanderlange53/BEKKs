test_that("Testing Asymmetric GridSearch works with StocksBonds", {
  set.seed(56285)
  signs <- as.matrix(c(-1,1))
  gs1 <- random_grid_search_asymmetric_BEKK(StocksBonds, 1 ,signs)

  # Starting values should be equal for same seed

  expect_equal(trunc(sum(gs1[[1]])), 3)
})


test_that("Testing Asymmetric GridSearch works with GoldStocksBonds", {
  set.seed(5463782)
  signs <- as.matrix(c(-1,1,-1))
  gs1 <- random_grid_search_asymmetric_BEKK(GoldStocksBonds, 1 ,signs)

  # Starting values should be equal for same seed

  expect_equal(trunc(sum(gs1[[1]])), 4)
})
