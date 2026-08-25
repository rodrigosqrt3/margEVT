# Helper: fit with one covariate
make_fit_with_cov <- function(seed = 1L) {
  set.seed(seed)
  n  <- 300L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.5), stats::runif(20, 5, 15)),
    mei  = stats::rnorm(n + 20L),
    tsa  = stats::rnorm(n + 20L),
    cos1 = cos(2 * pi * seq_len(n + 20L) / 365.25),
    sen1 = sin(2 * pi * seq_len(n + 20L) / 365.25),
    cos2 = cos(4 * pi * seq_len(n + 20L) / 365.25),
    sen2 = sin(4 * pi * seq_len(n + 20L) / 365.25)
  )
  fit_nhpp(df, threshold = 4,
           loc_vars  = c("mei", "tsa", "cos1", "sen1", "cos2", "sen2"),
           free_vars = c("cos1", "sen1", "cos2", "sen2"),
           lambda = 0, verbose = FALSE)
}

test_that("build_cov_annual returns correct number of rows", {
  fit <- make_fit_with_cov()
  df  <- build_cov_annual(fit, cov_vals = list(mei = 0, tsa = 0), n_obs = 365L)
  expect_equal(nrow(df), 365L)
})

test_that("seasonal columns are always present", {
  fit <- make_fit_with_cov()
  df  <- build_cov_annual(fit, cov_vals = list(mei = 0, tsa = 0))
  expect_true(all(c("cos1", "sen1", "cos2", "sen2") %in% names(df)))
})

test_that("build_cov_annual inherits temporal resolution from the fit", {
  fit <- make_fit_with_cov()
  fit$obs_per_year <- 52
  df <- build_cov_annual(fit, cov_vals = list(mei = 0, tsa = 0))
  expect_equal(nrow(df), 52L)
  expect_equal(df$cos1[1L], cos(2 * pi / 52))
})

test_that("scalar cov_vals are replicated to n_obs rows", {
  fit <- make_fit_with_cov()
  df  <- build_cov_annual(fit, cov_vals = list(mei = 1.5, tsa = -0.3),
                          n_obs = 365L)
  expect_true(all(df$mei == 1.5))
  expect_true(all(df$tsa == -0.3))
})

test_that("vector cov_vals of length n_obs are passed through", {
  fit <- make_fit_with_cov()
  v   <- stats::rnorm(365L)
  df  <- build_cov_annual(fit, cov_vals = list(mei = v, tsa = 0), n_obs = 365L)
  expect_equal(df$mei, v)
})

test_that("missing cov_vals fills with zero and warns", {
  fit <- make_fit_with_cov()
  expect_warning(
    df <- build_cov_annual(fit, cov_vals = list(), n_obs = 365L),
    regexp = "filled with 0"
  )
  expect_true(all(df$mei == 0))
})

test_that("interactions are computed correctly", {
  fit  <- make_fit_with_cov()
  ints <- list(mei_x_tsa = c("mei", "tsa"))
  df   <- build_cov_annual(fit,
                           cov_vals     = list(mei = 2, tsa = 3),
                           interactions = ints,
                           n_obs        = 10L)
  expect_true("mei_x_tsa" %in% names(df))
  expect_true(all(df$mei_x_tsa == 6))
})

test_that("interaction with missing column throws error", {
  fit  <- make_fit_with_cov()
  ints <- list(bad = c("mei", "nonexistent"))
  expect_error(
    suppressWarnings(
      build_cov_annual(fit, cov_vals = list(mei = 0), interactions = ints)
    ),
    regexp = "nonexistent"
  )
})

test_that("active_covariates returns NULL for stationary model", {
  set.seed(1L)
  df  <- data.frame(y = c(rexp(200, 0.5), runif(20, 5, 15)))
  fit <- fit_nhpp(df, threshold = 4, lambda = 0, verbose = FALSE)
  expect_null(active_covariates(fit))
})

test_that("active_covariates returns covariate names for non-stationary model", {
  fit <- make_fit_with_cov()
  ac  <- active_covariates(fit, tol = 1e-8)
  expect_true(all(c("mei", "tsa") %in% ac))
  expect_false(any(c("cos1", "sen1", "cos2", "sen2") %in% ac))
})

test_that("active_covariates inherits the fitted activity tolerance", {
  fit <- make_fit_with_cov()
  fit$active_tol <- 1e-2
  fit$par["mu.mei"] <- 5e-3
  fit$par["mu.tsa"] <- 2e-2
  ac <- active_covariates(fit)
  expect_false("mei" %in% ac)
  expect_true("tsa" %in% ac)
})

test_that("active_covariates supports legacy fits without activity tolerance", {
  fit <- make_fit_with_cov()
  fit$active_tol <- NULL
  fit$par["mu.mei"] <- 5e-3
  fit$par["mu.tsa"] <- 2e-2

  ac <- active_covariates(fit)

  expect_false("mei" %in% ac)
  expect_true("tsa" %in% ac)
})

test_that("active_covariates validates its activity tolerance", {
  fit <- make_fit_with_cov()
  expect_error(active_covariates(fit, tol = 0), regexp = "tol")
})

test_that("build_cov_annual and active_covariates reject non-nhpp_fit objects", {
  expect_error(build_cov_annual(list(), list()), regexp = "must be an nhpp_fit object")
  expect_error(active_covariates(list()), regexp = "must be an nhpp_fit object")
})

test_that("build_cov_annual validates interaction containers", {
  fit <- make_fit_with_cov()

  expect_error(
    build_cov_annual(fit, interactions = c("mei", "tsa")),
    regexp = "must be a named list"
  )
  expect_error(
    build_cov_annual(fit, interactions = list(c("mei", "tsa"))),
    regexp = "must be a named list"
  )

  duplicated <- list(a = c("mei", "tsa"), b = c("mei", "tsa"))
  names(duplicated) <- c("same", "same")
  expect_error(
    build_cov_annual(fit, interactions = duplicated),
    regexp = "unique names"
  )
})

test_that("build_cov_annual validates annual size and seasonal period", {
  fit <- make_fit_with_cov()
  expect_error(
    build_cov_annual(fit, n_obs = 0L),
    regexp = "n_obs.*positive integer"
  )
  expect_error(
    build_cov_annual(fit, period = 0),
    regexp = "period.*positive numeric"
  )
})

test_that("interaction inputs may use an existing seasonal column", {
  fit <- make_fit_with_cov()
  result <- build_cov_annual(
    fit,
    cov_vals = list(mei = 2, tsa = 0),
    n_obs = 10L,
    interactions = list(mei_x_cos1 = c("mei", "cos1"))
  )
  expect_equal(result$mei_x_cos1, 2 * result$cos1)
})

test_that("build_cov_annual rejects silently recyclable covariate lengths", {
  fit <- make_fit_with_cov()
  expect_error(
    build_cov_annual(fit, list(mei = 1:2, tsa = 0), n_obs = 10L),
    regexp = "scalar or a vector of length"
  )
})

test_that("build_cov_annual interaction validation throws error if not exactly 2 columns", {
  fit <- make_fit_with_cov()

  expect_error(
    suppressWarnings(
      build_cov_annual(fit, cov_vals = list(mei = 0), interactions = list(bad = c("mei")))
    ),
    regexp = "must name exactly 2 columns"
  )

  expect_error(
    suppressWarnings(
      build_cov_annual(fit, cov_vals = list(mei = 0), interactions = list(bad = c("mei", "tsa", "cos1")))
    ),
    regexp = "must name exactly 2 columns"
  )
})

test_that("active_covariates evaluates scale and shape parameters correctly", {
  set.seed(42L)
  n <- 300L
  df <- data.frame(
    y = c(stats::rexp(n, 0.5), stats::runif(20L, 5, 15)),
    scale_cov = stats::rnorm(n + 20L),
    shape_cov = stats::rnorm(n + 20L)
  )

  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = NULL,
                  scale_vars = "scale_cov",
                  shape_vars = "shape_cov",
                  penalty = "none", lambda = 0, verbose = FALSE)

  # Test block extraction deterministically rather than relying on random
  # covariates producing non-zero estimates in one optimizer run.
  fit$par["sigma.scale_cov"] <- 0.2
  fit$par["xi.shape_cov"] <- -0.1

  ac <- active_covariates(fit, tol = 1e-8)

  expect_true("scale_cov" %in% ac)
  expect_true("shape_cov" %in% ac)
})
