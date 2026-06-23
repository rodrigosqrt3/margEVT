# Helper: fit a simple model we can predict from
make_fit <- function(seed = 1L) {
  set.seed(seed)
  n  <- 300L
  df <- data.frame(
    y = c(stats::rexp(n, 0.5), stats::runif(20, 5, 15)),
    x = stats::rnorm(n + 20L),
    z = stats::rnorm(n + 20L)
  )
  fit_nhpp(df, threshold = 4, loc_vars = "x",
           lambda = 0, verbose = FALSE)
}

test_that("predict_params with newdata=NULL returns fitted values", {
  fit <- make_fit()
  p   <- predict_params(fit)
  expect_identical(p, fit$fitted)
})

test_that("predict_params returns correct length on newdata", {
  fit     <- make_fit()
  newdata <- data.frame(x = rnorm(10L))
  p       <- predict_params(fit, newdata)
  expect_length(p$mu,    10L)
  expect_length(p$sigma, 10L)
  expect_length(p$xi,    10L)
})

test_that("sigma is always positive", {
  fit     <- make_fit()
  newdata <- data.frame(x = rnorm(50L))
  p       <- predict_params(fit, newdata)
  expect_true(all(p$sigma > 0))
})

test_that("stationary model returns constant mu and sigma", {
  set.seed(1L)
  df  <- data.frame(y = c(rexp(200, 0.5), runif(20, 5, 15)))
  fit <- fit_nhpp(df, threshold = 4, lambda = 0, verbose = FALSE)
  p   <- predict_params(fit, data.frame(dummy = 1:5))
  expect_equal(length(unique(round(p$mu,    8L))), 1L)
  expect_equal(length(unique(round(p$sigma, 8L))), 1L)
  expect_equal(length(unique(round(p$xi,    8L))), 1L)
})

test_that("missing column in newdata throws informative error", {
  fit     <- make_fit()   # model uses column x
  newdata <- data.frame(z = rnorm(5L))   # x is missing
  expect_error(predict_params(fit, newdata), regexp = "not found in newdata")
})

test_that("predict_params rejects non nhpp_fit input", {
  expect_error(predict_params(list(a = 1)), regexp = "nhpp_fit")
})
