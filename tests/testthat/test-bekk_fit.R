# test_that("Symmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), -7382)
#   expect_equal(round(sum(x1$theta)),  2)
#   expect_equal(round(sum(c(x1$C0_t, x1$A_t, x1$G_t))), 1208)
# })
#
# test_that("Symmetric BEKK 2-dims works, ts object and QML_t_ratios = TRUE", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), -7382)
#   expect_equal(round(sum(x1$theta)),  2)
#   expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 1208)
# })
#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = FALSE", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), 75249)
#   expect_equal(round(sum(x1$theta)),  4)
#   expect_equal(round(sum(c(x1$C0_t, x1$A_t, x1$G_t))), 2384)
# })
#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = TRUE", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), 75249)
#   expect_equal(round(sum(x1$theta)),  4)
#   expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 2384)
# })
#
# test_that("Asymmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE", {
#   obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), -7340)
#   expect_equal(round(sum(x1$theta)),  3)
#   expect_equal(round(sum(c(x1$C0_t, x1$A_t, x1$G_t))), 1086)
# })
#
# test_that("Asymmetric BEKK 2-dims works, ts object and QML_t_ratios = TRUE", {
#   obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), -7340)
#   expect_equal(round(sum(x1$theta)),  3)
#   expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 1086)
# })
#
#
# test_that("Asymmetric BEKK 3-dims works, xts object and QML_t_ratios = TRUE", {
#   obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)
#
#   expect_equal(round(x1$log_likelihood), 75327)
#   expect_equal(round(sum(x1$theta)),  4)
#   expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 1797)
# })
test_that("Symmetric dBEKK 2-dims works, ts object and QML_t_ratios = FALSE", {
  obj_spec <- bekk_spec(model = list(type = "dbekk", asymmetric = F))
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  expect_equal(round(x1$log_likelihood), -7388)
  expect_equal(round(sum(x1$theta)),  3)
  expect_equal(round(sum(c(x1$C0_t, x1$A_t, x1$G_t))), 1233)
})
test_that("Asymmetric dBEKK 3-dims works, xts object and QML_t_ratios = TRUE", {
  obj_spec <- bekk_spec(model = list(type = "dbekk", asymmetric =T))
  x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)

  expect_equal(round(x1$log_likelihood), 75068)
  expect_equal(round(sum(x1$theta)),  4)
  expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 100)
})

test_that("Symmetric sBEKK 2-dims works, ts object and QML_t_ratios = FALSE", {
  obj_spec <- bekk_spec(model = list(type = "sbekk", asymmetric = F))
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  expect_equal(round(x1$log_likelihood), -7417)
  expect_equal(round(sum(x1$theta)),  1)
  expect_equal(round(sum(c(x1$C0_t, x1$A_t, x1$G_t))), 36)
})
test_that("Symmetric sBEKK 3-dims works, xts object and QML_t_ratios = TRUE", {
  obj_spec <- bekk_spec(model = list(type = "sbekk", asymmetric = F))
  x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = TRUE, max_iter = 50, crit = 1e-9)

  expect_equal(round(x1$log_likelihood), 75130)
  expect_equal(round(sum(x1$theta)),  4)
  expect_lt(sum(c(x1$C0_t, x1$A_t, x1$G_t)), 399.75)
})



