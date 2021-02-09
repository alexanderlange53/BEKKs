
test_that("Testing GridSearch works with TS_new", {
  set.seed(6723217)
  gs1 <- random_grid_search_BEKK(TS_Example,250)

  # Starting values should be equal for same seed
  set.seed(6723217)
  gs2 <- random_grid_search_BEKK(TS_Example,250)
  expect_identical(gs1[[1]], gs2[[1]])

})
