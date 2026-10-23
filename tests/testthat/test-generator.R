# =============================================================================
# test-generator.R
# Tests for fit_var_generator() and simulate_covariates() (R/generator.R)
# =============================================================================

library(testthat)
library(margEVT)

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

  base1 <- 2 * s$cos1 + 0.5 * s$sen2 + arima.sim(list(ar = 0.6), n)
  base2 <- 0.4 * base1 + s$sen1 + arima.sim(list(ar = 0.4), n)

  df$x1 <- as.numeric(base1)
  df$x2 <- as.numeric(base2)

  if (n_vars >= 3L) {
    df$x3 <- 3 * s$cos1 - 1.5 * s$sen2
  }

  df
}

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

test_that("fit_var_generator validates period and lag limits", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2"), period = 0),
    "period"
  )
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2"), lag_max = 0),
    "lag_max"
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

  testthat::local_mocked_bindings(active_covariates = function(...) NULL)

  expect_warning(
    result <- fit_var_generator(fit, df),
    "no active covariates"
  )
  expect_null(result)
})

test_that("fit_var_generator uses active_covariates() when vars = NULL", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()

  testthat::local_mocked_bindings(active_covariates = function(...) c("x1", "x2", "ghost_col"))

  gen <- fit_var_generator(fit, df)
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

test_that("fit_var_generator builds cos1/sen1/cos2/sen2 when absent from data", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(with_seasonal_cols = FALSE)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_s3_class(gen, "nhpp_var_generator")
  expect_true(all(c("cos1", "sen1", "cos2", "sen2") %in%
                    names(gen$seasonal_models$x1$coefs)))
})

test_that("fit_var_generator inherits the fitted temporal resolution", {
  fit <- .make_fake_fit()
  fit$obs_per_year <- 52
  df <- .make_test_data(n = 260L, with_seasonal_cols = FALSE)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_equal(gen$period, 52)
})

test_that("fit_var_generator reuses existing seasonal columns when present", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(with_seasonal_cols = TRUE)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_s3_class(gen, "nhpp_var_generator")
  expect_identical(gen$vars, c("x1", "x2"))
  expect_true(gen$spectral_radius < 1)
  expect_equal(gen$spectral_radius, max(gen$root_moduli))
})

test_that("fit_var_generator uses a genuine univariate AR model", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(n_vars = 2L)
  gen <- fit_var_generator(fit, df, vars = "x1")
  expect_s3_class(gen, "nhpp_var_generator")
  expect_identical(gen$vars, "x1")
  expect_identical(gen$K, 1L)
  expect_null(gen$fit_var)
  expect_equal(nrow(gen$A), 1L)
  expect_true(gen$spectral_radius < 1)
})

test_that("univariate generator fitting preserves the caller RNG stream", {
  fit <- .make_fake_fit()
  df <- .make_test_data(n_vars = 2L)
  set.seed(987L)
  expected <- runif(2L)

  set.seed(987L)
  first <- runif(1L)
  fit_var_generator(fit, df, vars = "x1")
  second <- runif(1L)

  expect_equal(c(first, second), expected)
})

test_that("companion matrix includes lag-shift rows for higher-order models", {
  A <- matrix(c(0.5, 0.2), nrow = 1L)
  result <- margEVT:::.companion_matrix(A, K = 1L, p = 2L)

  expect_equal(result, rbind(c(0.5, 0.2), c(1, 0)))
})

test_that("stationary covariance returns NULL for unsolvable or non-finite cases", {
  expect_null(
    margEVT:::.stationary_state_covariance(
      A = matrix(1, 1L, 1L), Sigma = matrix(1, 1L, 1L),
      K = 1L, p = 1L
    )
  )

  covariance_with_nonfinite_solution <-
    margEVT:::.stationary_state_covariance
  environment(covariance_with_nonfinite_solution) <- list2env(
    list(
      solve = function(lhs, rhs) rep(Inf, length(rhs)),
      .companion_matrix = margEVT:::.companion_matrix
    ),
    parent = environment(margEVT:::.stationary_state_covariance)
  )
  expect_null(
    covariance_with_nonfinite_solution(
      A = matrix(0, 1L, 1L), Sigma = matrix(1, 1L, 1L),
      K = 1L, p = 1L
    )
  )
})

test_that("fit_var_generator guards against near-zero residual sd", {
  fit <- .make_fake_fit()
  df  <- .make_test_data(n_vars = 3L)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x3"))
  expect_equal(gen$seasonal_models$x3$sd_res, 1)
})

test_that("fit_var_generator respects an explicit lag_max", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"), lag_max = 2L)
  expect_s3_class(gen, "nhpp_var_generator")
  expect_lte(gen$p_opt, 2L)
})

test_that("fit_var_generator rejects a non-finite selected lag", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  testthat::local_mocked_bindings(
    .select_var_order = function(...) NA_real_
  )
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2")),
    "did not return a finite order"
  )
})

test_that("fit_var_generator rejects an unstable selected VAR", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  testthat::local_mocked_bindings(
    .var_root_moduli = function(...) c(0.8, 1.01)
  )
  expect_error(
    fit_var_generator(fit, df, vars = c("x1", "x2")),
    "selected VAR.*is unstable"
  )
})

test_that("fit_var_generator computes a default lag_max when not supplied", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  expect_true(gen$p_opt >= 1L)
})

test_that("simulate_covariates returns NULL when generator is NULL", {
  expect_null(simulate_covariates(NULL, n_mc = 5L))
})

test_that("simulate_covariates validates the generator class", {
  expect_error(
    simulate_covariates(list(), n_mc = 5L),
    "must come from fit_var_generator"
  )
})

test_that("simulate_covariates validates simulation dimensions", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  expect_error(simulate_covariates(gen, n_mc = 0), "n_mc")
  expect_error(simulate_covariates(gen, n_mc = 1, burn_in = -1), "burn_in")
  expect_error(simulate_covariates(gen, n_mc = 1, n_obs = 0), "n_obs")
  expect_error(simulate_covariates(gen, n_mc = 1,
                                   initialization = "unknown"), "arg")
})

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

test_that("simulate_covariates inherits trajectory length from the generator", {
  fit <- .make_fake_fit()
  fit$obs_per_year <- 52
  df <- .make_test_data(n = 260L, with_seasonal_cols = FALSE)
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  mc <- simulate_covariates(gen, n_mc = 2L, burn_in = 10L, seed = 1L)
  expect_equal(vapply(mc, nrow, integer(1L)), c(52L, 52L))
})

test_that("simulate_covariates is reproducible with the same seed", {
  fit <- .make_fake_fit()
  df  <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  mc_a <- simulate_covariates(gen, n_mc = 2L, n_obs = 20L, burn_in = 10L, seed = 42L)
  mc_b <- simulate_covariates(gen, n_mc = 2L, n_obs = 20L, burn_in = 10L, seed = 42L)

  expect_equal(mc_a, mc_b)
})

test_that("simulate_covariates supports legacy generators without stored Sigma", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  gen$Sigma <- NULL

  mc <- simulate_covariates(
    gen, n_mc = 1L, n_obs = 10L,
    initialization = "burnin", burn_in = 5L, seed = 1L
  )
  expect_length(mc, 1L)
  expect_equal(nrow(mc[[1L]]), 10L)
})

test_that("simulate_covariates reports stationary initialization failure", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))

  testthat::local_mocked_bindings(
    .stationary_state_covariance = function(...) NULL
  )
  expect_error(
    simulate_covariates(
      gen, n_mc = 1L, n_obs = 10L, initialization = "stationary"
    ),
    regexp = "stationary initialization failed"
  )
})

test_that("stationary and burn-in initialization are both available", {
  fit <- .make_fake_fit()
  df <- .make_test_data()
  gen <- fit_var_generator(fit, df, vars = c("x1", "x2"))
  stationary <- simulate_covariates(
    gen, n_mc = 2L, n_obs = 20L,
    initialization = "stationary", seed = 1L
  )
  burnin <- simulate_covariates(
    gen, n_mc = 2L, n_obs = 20L, burn_in = 50L,
    initialization = "burnin", seed = 1L
  )
  expect_equal(vapply(stationary, nrow, integer(1L)), c(20L, 20L))
  expect_equal(vapply(burnin, nrow, integer(1L)), c(20L, 20L))
  expect_false(isTRUE(all.equal(stationary, burnin)))
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

# -----------------------------------------------------------------------------
# Least-squares VAR, deterministic and derived covariates, legacy objects
# -----------------------------------------------------------------------------

.make_var_data <- function(n = 600L, seed = 11L) {
  set.seed(seed)
  t <- seq_len(n)
  h <- seasonal_harmonics(t)
  x1 <- as.numeric(stats::arima.sim(list(ar = 0.6), n)) + 2 * h$cos1
  x2 <- 0.4 * x1 + as.numeric(stats::arima.sim(list(ar = 0.3), n)) + h$sen1
  cbind(data.frame(x1 = x1, x2 = x2), h)
}

test_that("least-squares VAR matches vars::VARselect and vars::VAR", {
  skip_if_not_installed("vars")
  df <- .make_var_data()
  gen <- fit_var_generator(.make_fake_fit(), df, vars = c("x1", "x2"), lag_max = 4L)
  Y <- sapply(c("x1", "x2"), function(v) {
    m <- gen$seasonal_models[[v]]
    res <- df[[v]] - as.numeric(cbind(1, as.matrix(df[gen$harmonics])) %*% m$coefs)
    (res - m$mu_res) / m$sd_res
  })
  sel <- vars::VARselect(Y, lag.max = 4L, type = "none")$selection[["SC(n)"]]
  expect_equal(gen$p_opt, as.integer(sel))
  ref <- vars::VAR(Y, p = gen$p_opt, type = "none")
  expect_equal(unname(gen$A), unname(vars::Bcoef(ref)), tolerance = 1e-8)
  expect_equal(sort(gen$root_moduli), sort(vars::roots(ref)), tolerance = 1e-8)
})

test_that("VAR helpers report rank deficiency and short samples", {
  Y <- cbind(a = stats::rnorm(50), b = 0)
  expect_true(is.na(margEVT:::.select_var_order(Y, 2L)))
  expect_null(margEVT:::.fit_var_ols(Y, 1L))
  expect_null(margEVT:::.fit_var_ols(matrix(stats::rnorm(6), 3), 2L))
  Yd <- cbind(a = stats::rnorm(50), b = stats::rnorm(50))
  Yd[, "b"] <- Yd[, "a"]
  expect_true(is.na(margEVT:::.select_var_order(Yd[-(1:2), ], 1L)))
})

test_that("an unidentifiable selected VAR is rejected", {
  df <- .make_var_data()
  testthat::local_mocked_bindings(.fit_var_ols = function(...) NULL)
  expect_error(fit_var_generator(.make_fake_fit(), df, vars = c("x1", "x2")),
               "not identifiable")
})

test_that("deterministic variables are excluded from the VAR and reconstructed", {
  df <- .make_var_data()
  df$season <- 3 * df$cos1 - df$sen2
  gen <- fit_var_generator(.make_fake_fit(), df, vars = c("x1", "season"))
  expect_identical(gen$stochastic_vars, "x1")
  expect_identical(gen$K, 1L)
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 30L, seed = 1L)
  expected <- 3 * cos(2 * pi * (1:30) / 365.25) - sin(4 * pi * (1:30) / 365.25)
  expect_equal(mc[[1L]]$season, expected, tolerance = 1e-8)
  expect_equal(mc[[2L]]$season, expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(mc[[1L]]$x1, mc[[2L]]$x1)))
})

test_that("a fully deterministic generator simulates without randomness", {
  df <- .make_var_data()
  df$season <- 2 + df$cos1
  gen <- fit_var_generator(.make_fake_fit(), df, vars = "season")
  expect_identical(gen$K, 0L)
  expect_identical(gen$p_opt, 0L)
  set.seed(5L); before <- stats::runif(1L)
  set.seed(5L)
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 10L)
  expect_identical(stats::runif(1L), before)
  expect_equal(mc[[1L]], mc[[2L]])
})

test_that("derive recomputes lagged covariates from simulated base variables", {
  df <- .make_var_data()
  derive <- function(b) data.frame(x1 = b$x1,
                                   x1_lag1 = c(NA, utils::head(b$x1, -1L)),
                                   x1_ma3 = as.numeric(stats::filter(
                                     b$x1, rep(1 / 3, 3), sides = 1)))
  gen <- fit_var_generator(.make_fake_fit(), df, vars = "x1", derive = derive,
                           derive_lead = 2L)
  mc <- simulate_covariates(gen, n_mc = 3L, n_obs = 40L, seed = 2L)
  d <- mc[[1L]]
  expect_identical(names(d), c("x1", "x1_lag1", "x1_ma3"))
  expect_equal(nrow(d), 40L)
  expect_equal(d$x1_lag1[-1L], d$x1[-40L])
  expect_equal(d$x1_ma3[3:40], (d$x1[1:38] + d$x1[2:39] + d$x1[3:40]) / 3)
})

test_that("derive output is validated", {
  df <- .make_var_data()
  gen_rows <- fit_var_generator(.make_fake_fit(), df, vars = "x1",
                                derive = function(b) b[-1L, , drop = FALSE])
  expect_error(simulate_covariates(gen_rows, n_mc = 1L, n_obs = 10L, seed = 1L),
               "must return a data frame with 10 rows")
  gen_na <- fit_var_generator(
    .make_fake_fit(), df, vars = "x1",
    derive = function(b) data.frame(l = c(NA, utils::head(b$x1, -1L))))
  expect_error(simulate_covariates(gen_na, n_mc = 1L, n_obs = 10L, seed = 1L),
               "increase `derive_lead`")
  expect_error(fit_var_generator(.make_fake_fit(), df, vars = "x1", derive = 1),
               "derive")
  expect_error(fit_var_generator(.make_fake_fit(), df, vars = "x1", derive_lead = -1),
               "derive_lead")
  expect_error(fit_var_generator(.make_fake_fit(), df, vars = "x1", n_harmonics = 0),
               "n_harmonics")
})

test_that("existing sin-named harmonics are reused and n_harmonics is respected", {
  df <- .make_var_data()
  df <- df[, c("x1", "x2")]
  df <- cbind(df, seasonal_harmonics(seq_len(nrow(df)), n_harmonics = 3L,
                                     sine = "sin"))
  gen <- fit_var_generator(.make_fake_fit(), df, vars = "x1", n_harmonics = 3L)
  expect_identical(gen$harmonics, c("cos1", "sin1", "cos2", "sin2",
                                    "cos3", "sin3"))
  expect_length(simulate_covariates(gen, n_mc = 1L, n_obs = 5L, seed = 1L), 1L)
})

test_that("legacy generators without harmonic and stochastic fields simulate", {
  df <- .make_var_data()
  gen <- fit_var_generator(.make_fake_fit(), df, vars = c("x1", "x2"))
  gen$harmonics <- NULL
  gen$stochastic_vars <- NULL
  gen$derive_lead <- NULL
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 8L, seed = 3L)
  expect_equal(dim(mc[[2L]]), c(8L, 2L))
})

test_that("simulation is vectorized, reproducible and seed-preserving", {
  df <- .make_var_data()
  gen <- fit_var_generator(.make_fake_fit(), df, vars = c("x1", "x2"))
  set.seed(9L); expected <- stats::runif(2L)
  set.seed(9L); first <- stats::runif(1L)
  a <- simulate_covariates(gen, n_mc = 200L, n_obs = 365L, seed = 4L)
  second <- stats::runif(1L)
  expect_equal(c(first, second), expected)
  b <- simulate_covariates(gen, n_mc = 200L, n_obs = 365L, seed = 4L)
  expect_identical(a, b)
  expect_error(simulate_covariates(gen, n_mc = 1L, seed = "x"), "seed")
  # stationary marginal variance of the simulated anomalies matches the VAR
  x <- vapply(a, function(d) d$x1[365L], numeric(1L))
  expect_gt(stats::sd(x), 0)
})
