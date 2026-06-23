make_select_df <- function(n = 400L, seed = 7L) {
  set.seed(seed)
  data.frame(
    y = c(stats::rexp(n, 0.3), stats::runif(40L, 5, 20)),
    x = stats::rnorm(n + 40L),
    z = stats::rnorm(n + 40L)
  )
}

test_that("lambda=bic returns a converged fit", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = c("x", "z"),
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_s3_class(fit, "nhpp_fit")
  expect_true(fit$converged)
  expect_true(is.numeric(fit$lambda))
  expect_true(all(fit$lambda > 0))
})

test_that("bic selection returns a named lambda vector", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = "x",
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_named(fit$lambda, c("mu", "sigma", "xi"))
})

test_that("bic lasso shrinks noise covariate more than signal", {
  set.seed(99L)
  n   <- 500L
  # x is correlated with extremes, z is pure noise
  y   <- c(stats::rexp(n, 0.3), stats::runif(50L, 5, 20))
  x   <- c(stats::rnorm(n), stats::rnorm(50L, mean = 2))
  z   <- stats::rnorm(n + 50L)
  df  <- data.frame(y = y, x = x, z = z)
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = c("x", "z"),
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  # Signal covariate should have larger absolute coefficient than noise
  expect_gte(abs(fit$par["mu.x"]), abs(fit$par["mu.z"]))
})

test_that("ridge bic fit has non-zero coefficients for all covariates", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = c("x", "z"),
                  penalty  = "ridge",
                  lambda   = "bic",
                  verbose  = FALSE)
  # Ridge never shrinks exactly to zero
  expect_gt(abs(fit$par["mu.x"]), 1e-6)
  expect_gt(abs(fit$par["mu.z"]), 1e-6)
})
