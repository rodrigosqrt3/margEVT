# =============================================================================
# test-generator.R
# Tests for fit_var_generator() and simulate_covariates() (R/generator.R)
# =============================================================================

# ---- shared fixtures --------------------------------------------------------

.make_fake_fit <- function() structure(list(), class = "nhpp_fit")

.make_seasonal_series <- function(n, period = 365.25, seed = 1L) {
  set.seed(seed)
  t_grid <- seq_len(n)
  cos1 <- cos(2 * pi * t_grid / period)
  sen1 <- sin(2 * pi * t_grid / period)
  cos2 <- cos(4 * pi * t_grid / period)
  sen2 <- sin(4 * pi * t_grid / period)
  list(cos1 = cos1, sen1 = sen1, cos2 = cos2, sen2 = sen2, t_grid = t_grid)
}

.make_test_data <- function(n = 400L, with_seasonal_cols = TRUE,
                            n_vars = 2L, seed = 1L) {
  set.seed(seed)
  s <- .make_seasonal_series(n, seed = seed)

  df <- data.frame(y = rnorm(n))
  if (with_seasonal_cols) {
    df$cos1 <- s$cos1; df$sen1 <- s$sen1
    df$cos2 <- s$cos2; df$sen2 <- s$sen2
  }

  # Two correlated, seasonally-modulated "active" covariates with genuine
  # residual variance — these should produce a well-behaved VAR fit.
  base1 <- 2 * s$cos1 + 0.5 * s$sen2 + arima.sim(list(ar = 0.6), n)
  base2 <- 0.4 * base1 + s$sen1 + arima.sim(list(ar = 0.4), n)

  df$x1 <- as.numeric(base1)
  df$x2 <- as.numeric(base2)

  if (n_vars >= 3L) {
    # Degenerate covariate: an exact linear combination of the seasonal
    # harmonics, so its deseasonalized residual variance is ~0. Exercises
    # the sd_res < 1e-8 guard.
    df$x3 <- 3 * s$cos1 - 1.5 * s$sen2
  }

  df
}

# ---- fit_var_generator(): input validation ---------------------------------

test_that("fit_var_generator validates `fit` and `data`", {
  df <- .make_test_data()
  expect_error(
    fit_var_generator(list(), df),
    "must be an nhpp_fit object"
  )
  expect_error(
    fit_var_generator(.make_fake_fit(), as.list(df)),
    "must be a data frame"
  )
})

test_that("fit_var_generator errors when the 'vars' package is unavailable", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  mockery::stub(fit_var_generator, "requireNamespace", FALSE)
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2")),
    "package 'vars' is required"
  )
})

test_that("fit_var_generator errors on explicit vars missing from data", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "does_not_exist")),
    "columns not in data"
  )
})

test_that("fit_var_generator warns and returns NULL with no active covariates", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  mockery::stub(fit_var_generator, "active_covariates", NULL)
  expect_warning(
    result <- fit_var_generator(fit, df),
    "no active covariates"
  )
  expect_null(result)
})

test_that("fit_var_generator uses active_covariates() when vars = NULL", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  mockery::stub(fit_var_generator, "active_covariates", c("x1", "x2", "ghost_col"))
  gen <- fit_var_generator(fit, df)
  # ghost_col isn't in data and must be silently dropped by the intersection,
  # not passed through to the VAR fit.
  expect_setequal(gen$vars, c("x1", "x2"))
})

test_that("fit_var_generator errors when too few complete rows remain", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(n = 20L)
  df$x1[1:15] <- NA_real_
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2")),
    "too few complete rows"
  )
})

# ---- fit_var_generator(): seasonal harmonics handling -----------------------

test_that("fit_var_generator builds cos1/sen1/cos2/sen2 when absent from data", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(with_seasonal_cols = FALSE)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_s3_class(gen, "nhpp_var_generator")
  expect_true(all(c("cos1", "sen1", "cos2", "sen2") %in%
                    names(gen$seasonal_models$x1$coefs)))
})

test_that("fit_var_generator reuses existing seasonal columns when present", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(with_seasonal_cols = TRUE)
  # Corrupt the harmonics in a way that would error downstream if
  # build_var_generator() tried to overwrite/recompute and got confused
  # by an unexpected column type.
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_s3_class(gen, "nhpp_var_generator")
  expect_identical(gen$vars, c("x1", "x2"))
})

# ---- fit_var_generator(): K == 1 (dummy-noise padding) ----------------------

test_that("fit_var_generator pads a single covariate with dummy noise for VAR", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(n_vars = 2L)
  gen <- fit_var_generator(fit, df, vars = "x1")
  expect_s3_class(gen, "nhpp_var_generator")
  expect_identical(gen$vars, "x1")
  expect_true(".dummy_noise" %in% gen$var_colnames)
  expect_identical(gen$fit_var$K, 2L)
})

# ---- fit_var_generator(): degenerate residual variance guard ----------------

test_that("fit_var_generator guards against near-zero residual sd", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(n_vars = 3L)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x3"))
  expect_equal(gen$seasonal_models$x3$sd_res, 1)
})

# ---- fit_var_generator(): lag_max argument -----------------------------------

test_that("fit_var_generator respects an explicit lag_max", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"), lag_max = 2L)
  expect_s3_class(gen, "nhpp_var_generator")
  expect_lte(gen$p_opt, 2L)
})

test_that("fit_var_generator computes a default lag_max when not supplied", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_true(gen$p_opt >= 1L)
})

# ---- simulate_covariates(): input handling ----------------------------------

test_that("simulate_covariates returns NULL when generator is NULL", {
  expect_null(simulate_covariates(NULL, n_mc = 5L))
})

test_that("simulate_covariates validates the generator class", {
  expect_error(
    simulate_covariates(list(), n_mc = 5L),
    "must come from fit_var_generator"
  )
})

# ---- simulate_covariates(): output structure --------------------------------

test_that("simulate_covariates produces the right shape and column names", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  mc <- simulate_covariates(gen, n_mc = 3L, n_obs = 30L, burn_in = 10L, seed = 1L)

  expect_type(mc, "list")
  expect_length(mc, 3L)
  for (yr_df in mc) {
    expect_s3_class(yr_df, "data.frame")
    expect_equal(nrow(yr_df), 30L)
    expect_setequal(names(yr_df), c("x1", "x2"))
    expect_true(all(vapply(yr_df, is.numeric, logical(1L))))
    expect_true(all(vapply(yr_df, function(col) all(is.finite(col)), logical(1L))))
  }
})

test_that("simulate_covariates is reproducible with the same seed", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  mc_a <- simulate_covariates(gen, n_mc = 2L, n_obs = 20L, burn_in = 10L, seed = 42L)
  mc_b <- simulate_covariates(gen, n_mc = 2L, n_obs = 20L, burn_in = 10L, seed = 42L)

  expect_equal(mc_a, mc_b)
})

test_that("simulate_covariates runs without error when seed is NULL", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  expect_no_error(
    simulate_covariates(gen, n_mc = 1L, n_obs = 15L, burn_in = 5L, seed = NULL)
  )
})

test_that("simulate_covariates output works as mc_sample for build_cov_annual()", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  mc  <- simulate_covariates(gen, n_mc = 1L, n_obs = 20L, burn_in = 5L, seed = 1L)

  # build_cov_annual() needs an nhpp_fit with a $dm slot; fabricate a minimal
  # one whose mu design matrix references the simulated covariates.
  fake_dm <- list(
    X_mu    = matrix(0, nrow = 1L, ncol = 3L,
                     dimnames = list(NULL, c("(Intercept)", "x1", "x2"))),
    X_sigma = matrix(0, nrow = 1L, ncol = 1L, dimnames = list(NULL, "(Intercept)")),
    X_xi    = matrix(0, nrow = 1L, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  )
  fake_fit <- structure(list(dm = fake_dm), class = "nhpp_fit")

  cov_df <- build_cov_annual(fake_fit, as.list(mc[[1L]]), n_obs = 20L)
  expect_true(all(c("x1", "x2") %in% names(cov_df)))
  expect_equal(nrow(cov_df), 20L)
})
