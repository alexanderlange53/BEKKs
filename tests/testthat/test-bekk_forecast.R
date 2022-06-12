test_that("Symmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE, n.ahead = 1", {
  obj_spec <- bekk_spec()
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_forecast(x1, n.ahead = 1)

  expect_equal(round(as.numeric(x2$volatility_forecast[1]), 2), 0.2)
  expect_equal(round(as.numeric(x2$volatility_forecast[2]), 2), -0.2)
  expect_equal(round(as.numeric(x2$volatility_forecast[3]), 2), 0.6)

  expect_equal(round(sum(x2$H_t_forecast), 2), 0.35)
  expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 0.59)
  expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 0.62)
})

test_that("Symmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE, n.ahead = 5", {
  obj_spec <- bekk_spec()
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_forecast(x1, n.ahead = 5)

  expect_equal(round(sum(x2$volatility_forecast[,1]), 2), 1)
  expect_equal(round(sum(x2$volatility_forecast[,2]), 2), -0.99)
  expect_equal(round(sum(x2$volatility_forecast[,3]), 3), 3.055)

  expect_equal(round(sum(x2$H_t_forecast), 2), 1.83)
  expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 2.98)
  expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 3.17)
})

#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = FALSE, n.ahead = 1", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   x2 <- bekk_forecast(x1, n.ahead = 1)
#
#   expect_equal(round(as.numeric(x2$volatility_forecast[,1]), 2), 0.01)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,2]), 2), 0.12)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,3]), 2), 0.22)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,4]), 2), 0.01)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,5]), 2), -0.16)
#
#   expect_equal(round(sum(x2$H_t_forecast), 2), 0)
#   expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 0.19)
#   expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 0.23)
# })
#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = FALSE, n.ahead = 5", {
#   obj_spec <- bekk_spec()
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   x2 <- bekk_forecast(x1, n.ahead = 5)
#
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,1])), 2), 0.05)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,2])), 2), 0.58)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,3])), 2), 1.09)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,4])), 2), 0.04)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,5])), 2), -0.8)
#
#   expect_equal(round(sum(x2$H_t_forecast), 2), 0)
#   expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 0.85)
#   expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 1.15)
# })


test_that("Asymmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE, n.ahead = 1", {
  obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_forecast(x1, n.ahead = 1)

  expect_equal(round(as.numeric(x2$volatility_forecast[1]), 2), 0.2)
  expect_equal(round(as.numeric(x2$volatility_forecast[2]), 2), -0.24)
  expect_equal(round(as.numeric(x2$volatility_forecast[3]), 2), 0.59)

  expect_equal(round(sum(x2$H_t_forecast), 2), 0.33)
  expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 0.51)
  expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 0.57)
})

test_that("Symmetric BEKK 2-dims works, ts object and QML_t_ratios = FALSE, n.ahead = 5", {
  obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_forecast(x1, n.ahead = 5)

  expect_equal(round(sum(as.numeric(x2$volatility_forecast[,1])), 2), 1)
  expect_equal(round(sum(as.numeric(x2$volatility_forecast[,2])), 2), -1.16)
  expect_equal(round(sum(as.numeric(x2$volatility_forecast[,3])), 2), 2.85)

  expect_equal(round(sum(x2$H_t_forecast), 2), 1.56)
  expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 2.57)
  expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 2.86)
})

#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = FALSE, n.ahead = 1", {
#   obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   x2 <- bekk_forecast(x1, n.ahead = 1)
#
#   expect_equal(round(as.numeric(x2$volatility_forecast[,1]), 2), 0.01)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,2]), 2), 0.01)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,3]), 2), 0.35)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,4]), 2), 0.01)
#   expect_equal(round(as.numeric(x2$volatility_forecast[,5]), 2), -0.11)
#
#   expect_equal(round(sum(x2$H_t_forecast), 2), 0)
#   expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 0.25)
#   expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 0.31)
# })
#
# test_that("Symmetric BEKK 3-dims works, xts object and QML_t_ratios = FALSE, n.ahead = 5", {
#   obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   x2 <- bekk_forecast(x1, n.ahead = 5)
#
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,1])), 2), 0.05)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,2])), 2), 0.09)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,3])), 2), 1.64)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,4])), 2), 0.05)
#   expect_equal(round(sum(as.numeric(x2$volatility_forecast[,5])), 2), -0.63)
#
#   expect_equal(round(sum(x2$H_t_forecast), 2), 0)
#   expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 1.04)
#   expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 1.45)
# })



test_that("Symmetric scalar BEKK 2-dims works, ts object and QML_t_ratios = FALSE, n.ahead = 1", {
  obj_spec <- bekk_spec(model = list(type="sbekk", asymmetric =F))
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_forecast(x1, n.ahead = 10)

  expect_equal(round(as.numeric(x2$volatility_forecast[1,1]), 2), 0.19)
  expect_equal(round(as.numeric(x2$volatility_forecast[2,2]), 2), -0.1)
  expect_equal(round(as.numeric(x2$volatility_forecast[2,3]), 2), 0.64)

  expect_equal(round(sum(x2$H_t_forecast), 2), 4.28)
  expect_equal(round(sum(x2$volatility_lower_conf_band), 2), 7.38)
  expect_equal(round(sum(x2$volatility_upper_conf_band), 2), 7.49)
})
