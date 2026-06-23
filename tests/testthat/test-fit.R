make_test_df <- function(n = 300L, seed = 42L) {
  set.seed(seed)
  data.frame(
    y = c(stats::rexp(n, rate = 0.5), stats::runif(30L, 5, 15)),
    x = stats::rnorm(n + 30L)
  )
}

test_that("fit_nhpp returns an nhpp_fit object", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_s3_class(fit, "nhpp_fit")
})

test_that("penalty=none gives lambda=0 and correct parameter names", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$lambda, 0)
  expect_named(fit$par,
               c("mu.(Intercept)", "sigma.(Intercept)", "xi.(Intercept)"))
})

test_that("penalty=lasso forces alpha=1", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 1)
  expect_equal(fit$penalty, "lasso")
})

test_that("penalty=ridge forces alpha=0", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "ridge", lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 0)
  expect_equal(fit$penalty, "ridge")
})

test_that("penalty=elnet uses supplied alpha", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "elnet", alpha = 0.3,
                  lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 0.3)
})

test_that("large lambda shrinks covariate to zero under lasso", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 500,
                  verbose = FALSE)
  expect_lt(abs(fit$par["mu.x"]), 0.01)
})

test_that("free_vars not shrunk by large lambda", {
  df       <- make_test_df()
  df$trend <- seq(0, 1, length.out = nrow(df))
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars  = c("trend", "x"),
                  free_vars = "trend",
                  penalty   = "lasso", lambda = 500,
                  verbose   = FALSE)
  expect_lt(abs(fit$par["mu.x"]),    0.01)
  expect_gt(abs(fit$par["mu.trend"]), 0)
})

test_that("numeric lambda skips grid search", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 0.5,
                  verbose = FALSE)
  expect_equal(fit$lambda, 0.5)
})

test_that("coef() returns parameter vector", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_identical(coef(fit), fit$par)
})

test_that("fitted values have correct length and sigma positive", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_length(fit$fitted$mu,    nrow(df))
  expect_length(fit$fitted$sigma, nrow(df))
  expect_length(fit$fitted$xi,    nrow(df))
  expect_true(all(fit$fitted$sigma > 0))
})

test_that("calc_hessian stores a square matrix", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none",
                  calc_hessian = TRUE, verbose = FALSE)
  expect_true(is.matrix(fit$hessian))
  expect_equal(nrow(fit$hessian), length(fit$par))
})
