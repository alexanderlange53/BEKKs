test_that("Testing Asymmetric GridSearch works with StocksBonds", {
  set.seed(56285)
  signs <- as.matrix(c(-1,1))
  gs1 <- random_grid_search_asymmetric_BEKK(StocksBonds, 1 ,signs)

  # Starting values should be equal for same seed

  expect_equal(sum(gs1[[1]]), sum(c(0.0005818369, -0.0014075557,  0.0040875551,  0.3036117923, -0.0001716152, -0.0002938880,  0.3003445156,  0.9358654946,  0.0004830816, -0.0002328148,  0.9310227285)))
})


test_that("Testing Asymmetric GridSearch works with GoldStocksBonds", {
  set.seed(5463782)
  signs <- as.matrix(c(-1,1,-1))
  gs1 <- random_grid_search_asymmetric_BEKK(GoldStocksBonds, 1 ,signs)

  # Starting values should be equal for same seed
  sum(gs1[[1]])
  expect_equal(sum(gs1[[1]]), 4.592864)
})
