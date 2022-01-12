test_that("Simulating 2 dim symmetric BEKK based on estimated ts object", {
  obj_spec <- bekk_spec()
  x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_sim(x1, nobs = 100)
  expect_equal(nrow(x2), 100)
  expect_equal(ncol(x2), 2)
})

test_that("Simulating 2 dim symmetric BEKK based on bekk_spec", {
  obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = FALSE), N = 2,
                        init_values = c(0.019783677, -0.014962718, 0.085265736, 0.178859630, -0.007691659,
                                         -0.037152006, 0.298147337, 0.981184605, 0.001337367, 0.013260377, 0.951431611))

  x2 <- bekk_sim(obj_spec, nobs = 100)
  expect_equal(nrow(x2), 100)
  expect_equal(ncol(x2), 2)
})

test_that("Simulating 3 dim symmetric BEKK based on estimated xts object", {
  obj_spec <- bekk_spec()
  x1 <- bekk_fit(obj_spec, GoldStocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)

  x2 <- bekk_sim(x1, nobs = 100)
  expect_equal(nrow(x2), 100)
  expect_equal(ncol(x2), 3)
})

test_that("Simulating 3 dim symmetric BEKK based on bekk_spec", {
  obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = FALSE), N = 3,
                        init_values = c(0.0010122220, -0.0002678085, 0.0001905036, 0.0008937830, 0.0003111280,
                                        0.0003893194, 0.2453276126, 0.0222127631, -0.0434456548, -0.0366629406,
                                        0.2637613982, -0.0248105673, 0.0624407934, -0.0016446789, 0.1733957064,
                                        0.9651528236, -0.0060973318, 0.0155928988, 0.0125157899, 0.9606411455,
                                        0.0008724460, -0.0196494968, -0.0009070321, 0.9801993180))

  x2 <- bekk_sim(obj_spec, nobs = 100)
  expect_equal(nrow(x2), 100)
  expect_equal(ncol(x2), 3)
})

# test_that("Simulating 2 dim asymmetric BEKK based on estimated ts object", {
#   obj_spec <- bekk_spec(model = list(type = 'bekk', asymmetric = TRUE))
#   x1 <- bekk_fit(obj_spec, StocksBonds, QML_t_ratios = FALSE, max_iter = 50, crit = 1e-9)
#
#   x2 <- bekk_sim(x1, nobs = 100)
#   expect_equal(nrow(x2), 100)
#   expect_equal(ncol(x2), 2)
# })

test_that("Simulating 2 dim asymmetric BEKK based on bekk_spec", {
  obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = TRUE), N = 2,
                        init_values = as.matrix(c(0.019783677, -0.014962718, 0.085265736, 0.178859630, -0.007691659,
                                        -0.037152006, 0.298147337, 0.05234,-0.02,0.032,0.1032, 0.981184605, 0.001337367, 0.013260377, 0.951431611)))

  x2 <- bekk_sim(obj_spec, nobs = 100)
  expect_equal(nrow(x2), 100)
  expect_equal(ncol(x2), 2)
})
