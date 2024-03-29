test_that("virf works with 2-dimensional test data set", {
  spec = bekk_spec()
  x = bekk_fit(spec, StocksBonds)
  s = virf(x)
  expect_equal(round(sum(s$VIRF),3), -0.727)
  expect_equal(round(sum(s$VIRF_upper),3), -0.661)
  expect_equal(round(sum(s$VIRF_lower),3), -0.794)
})
# test_that("virf works with 3-dimensional test data set", {
#   spec = bekk_spec()
#   x = bekk_fit(spec, GoldStocksBonds)
#   s = virf(x, n.ahead = 50)
#   expect_equal(round(sum(s$VIRF_upper),3), 0)
#   expect_equal(round(sum(s$VIRF_lower),3), 0)
#   expect_equal(round(sum(s$VIRF),3), 0)
# })
