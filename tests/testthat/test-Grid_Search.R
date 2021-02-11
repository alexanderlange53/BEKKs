test_that("Testing GridSearch works with TS_new", {
  set.seed(6723217)
  gs1 <- random_grid_search_BEKK(data.matrix(BI),250)

  # Starting values should be equal for same seed

  expect_equal(gs1[[1]], c(0.19663233, -0.22513552,  0.57296538,  0.03252316,  0.01099530,  0.53818486, -0.42037775, 0.59263607,-0.26044763,  0.98807132, -0.05697966))
})
