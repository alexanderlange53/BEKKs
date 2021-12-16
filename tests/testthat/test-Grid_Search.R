test_that("Testing GridSearch works with StocksBonds", {
  set.seed(6723217)
  gs1 <- random_grid_search_BEKK(data.matrix(StocksBonds))

  # Starting values should be equal for same seed

  expect_equal(sum(gs1[[1]]), sum(c(0.0005818369, -0.0014075557,  0.0040875551,  0.3036117923, -0.0001716152, -0.0002938880,  0.3003445156,  0.9358654946,  0.0004830816, -0.0002328148,  0.9310227285)))
})
