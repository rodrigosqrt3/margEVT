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

test_that("fit_nhpp input validation throws expected errors", {
  df <- make_test_df()
  expect_error(fit_nhpp(list(y = 1:10), threshold = 4), regexp = "must be a data frame")
  expect_error(fit_nhpp(data.frame(x = 1:10), threshold = 4), regexp = "must contain a column named `y`")
  expect_error(fit_nhpp(df, threshold = c(1, 4)), regexp = "single numeric value")
  expect_error(fit_nhpp(df, threshold = 4, penalty = "lasso", lambda = "invalid"), regexp = "positive numeric value or")
})

test_that("fit_nhpp throws warning on < 5 exceedances", {
  df_few <- data.frame(y = c(1, 2, 3, 4.5, 4.5, 4.5), x = rnorm(6))
  expect_warning(
    fit_nhpp(df_few, threshold = 4, penalty = "none", verbose = TRUE),
    regexp = "fewer than 5 exceedances"
  )
})

test_that("fit_nhpp throws warning on non-convergence", {
  df <- make_test_df()
  expect_warning(
    fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none", maxit = 1L, verbose = TRUE),
    regexp = "optimizer did not converge"
  )
})

# -----------------------------------------------------------------------------
# S3 Print Method Test
# -----------------------------------------------------------------------------
test_that("print.nhpp_fit prints without error", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)

  expect_output(print(fit), regexp = "-- nhpp_fit --")
  expect_output(print(fit), regexp = "Threshold")
  expect_output(print(fit), regexp = "Active params")
  expect_output(print(fit), regexp = "Coefficients")
})

test_that(".warmstart handles empty pen_idx_all", {
  df <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = NULL, scale_vars = NULL,
                  penalize_shape = FALSE, penalty = "lasso", lambda = 1.0, verbose = FALSE)
  expect_s3_class(fit, "nhpp_fit")
})

test_that(".fit_at_lambda handles missing init correctly", {
  df <- make_test_df()
  dm <- margEVT:::build_design_matrices(df, loc_vars = "x", scale_vars = NULL, shape_vars = NULL, free_vars = NULL)
  attr(dm, "threshold") <- 4

  res <- margEVT:::.fit_at_lambda(dm, df$y, threshold = 4, lambda = 0, alpha = 1,
                                  penalize_shape = TRUE, init = NULL)

  expect_type(res, "list")
  expect_true(res$converged)
})

test_that(".warmstart handles all(lambda == 0) early exit", {
  init_mock <- c(1.0, 2.0, 3.0)
  res <- margEVT:::.warmstart(init = init_mock, dm = NULL, y = NULL,
                              threshold = 0, lambda = 0, alpha = 1,
                              penalize_shape = TRUE)
  expect_identical(res, init_mock)
})
