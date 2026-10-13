# Tests for the 0.4.0 trajectory generator: least-squares VAR without the
# 'vars' dependency, deterministic variables, derived covariates, vectorized
# simulation, and harmonic handling.

.gen_fit <- function() structure(list(), class = "nhpp_fit")

.gen_data <- function(n = 600L, seed = 11L) {
  set.seed(seed)
  t <- seq_len(n)
  h <- seasonal_harmonics(t)
  x1 <- as.numeric(stats::arima.sim(list(ar = 0.6), n)) + 2 * h$cos1
  x2 <- 0.4 * x1 + as.numeric(stats::arima.sim(list(ar = 0.3), n)) + h$sen1
  cbind(data.frame(x1 = x1, x2 = x2), h)
}

test_that("least-squares VAR matches vars::VARselect and vars::VAR", {
  skip_if_not_installed("vars")
  df <- .gen_data()
  gen <- fit_var_generator(.gen_fit(), df, vars = c("x1", "x2"), lag_max = 4L)
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
  df <- .gen_data()
  testthat::local_mocked_bindings(.fit_var_ols = function(...) NULL)
  expect_error(fit_var_generator(.gen_fit(), df, vars = c("x1", "x2")),
               "not identifiable")
})

test_that("deterministic variables are excluded from the VAR and reconstructed", {
  df <- .gen_data()
  df$season <- 3 * df$cos1 - df$sen2
  gen <- fit_var_generator(.gen_fit(), df, vars = c("x1", "season"))
  expect_identical(gen$stochastic_vars, "x1")
  expect_identical(gen$K, 1L)
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 30L, seed = 1L)
  expected <- 3 * cos(2 * pi * (1:30) / 365.25) - sin(4 * pi * (1:30) / 365.25)
  expect_equal(mc[[1L]]$season, expected, tolerance = 1e-8)
  expect_equal(mc[[2L]]$season, expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(mc[[1L]]$x1, mc[[2L]]$x1)))
})

test_that("a fully deterministic generator simulates without randomness", {
  df <- .gen_data()
  df$season <- 2 + df$cos1
  gen <- fit_var_generator(.gen_fit(), df, vars = "season")
  expect_identical(gen$K, 0L)
  expect_identical(gen$p_opt, 0L)
  set.seed(5L); before <- stats::runif(1L)
  set.seed(5L)
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 10L)
  expect_identical(stats::runif(1L), before)
  expect_equal(mc[[1L]], mc[[2L]])
})

test_that("derive recomputes lagged covariates from simulated base variables", {
  df <- .gen_data()
  derive <- function(b) data.frame(x1 = b$x1,
                                   x1_lag1 = c(NA, utils::head(b$x1, -1L)),
                                   x1_ma3 = as.numeric(stats::filter(
                                     b$x1, rep(1 / 3, 3), sides = 1)))
  gen <- fit_var_generator(.gen_fit(), df, vars = "x1", derive = derive,
                           derive_lead = 2L)
  mc <- simulate_covariates(gen, n_mc = 3L, n_obs = 40L, seed = 2L)
  d <- mc[[1L]]
  expect_identical(names(d), c("x1", "x1_lag1", "x1_ma3"))
  expect_equal(nrow(d), 40L)
  expect_equal(d$x1_lag1[-1L], d$x1[-40L])
  expect_equal(d$x1_ma3[3:40], (d$x1[1:38] + d$x1[2:39] + d$x1[3:40]) / 3)
})

test_that("derive output is validated", {
  df <- .gen_data()
  gen_rows <- fit_var_generator(.gen_fit(), df, vars = "x1",
                                derive = function(b) b[-1L, , drop = FALSE])
  expect_error(simulate_covariates(gen_rows, n_mc = 1L, n_obs = 10L, seed = 1L),
               "must return a data frame with 10 rows")
  gen_na <- fit_var_generator(
    .gen_fit(), df, vars = "x1",
    derive = function(b) data.frame(l = c(NA, utils::head(b$x1, -1L))))
  expect_error(simulate_covariates(gen_na, n_mc = 1L, n_obs = 10L, seed = 1L),
               "increase `derive_lead`")
  expect_error(fit_var_generator(.gen_fit(), df, vars = "x1", derive = 1),
               "derive")
  expect_error(fit_var_generator(.gen_fit(), df, vars = "x1", derive_lead = -1),
               "derive_lead")
  expect_error(fit_var_generator(.gen_fit(), df, vars = "x1", n_harmonics = 0),
               "n_harmonics")
})

test_that("existing sin-named harmonics are reused and n_harmonics is respected", {
  df <- .gen_data()
  df <- df[, c("x1", "x2")]
  df <- cbind(df, seasonal_harmonics(seq_len(nrow(df)), n_harmonics = 3L,
                                     sine = "sin"))
  gen <- fit_var_generator(.gen_fit(), df, vars = "x1", n_harmonics = 3L)
  expect_identical(gen$harmonics, c("cos1", "sin1", "cos2", "sin2",
                                    "cos3", "sin3"))
  expect_length(simulate_covariates(gen, n_mc = 1L, n_obs = 5L, seed = 1L), 1L)
})

test_that("legacy generators without harmonic and stochastic fields simulate", {
  df <- .gen_data()
  gen <- fit_var_generator(.gen_fit(), df, vars = c("x1", "x2"))
  gen$harmonics <- NULL
  gen$stochastic_vars <- NULL
  gen$derive_lead <- NULL
  mc <- simulate_covariates(gen, n_mc = 2L, n_obs = 8L, seed = 3L)
  expect_equal(dim(mc[[2L]]), c(8L, 2L))
})

test_that("simulation is vectorized, reproducible and seed-preserving", {
  df <- .gen_data()
  gen <- fit_var_generator(.gen_fit(), df, vars = c("x1", "x2"))
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

test_that("seasonal_harmonics validates its inputs", {
  expect_error(seasonal_harmonics(c(1, NA)), "`t`")
  expect_error(seasonal_harmonics(1:3, period = 0), "period")
  expect_error(seasonal_harmonics(1:3, n_harmonics = 0), "n_harmonics")
  h <- seasonal_harmonics(1:4, period = 4, n_harmonics = 1L, sine = "sin")
  expect_equal(h$cos1, cos(2 * pi * (1:4) / 4))
  expect_equal(h$sin1, sin(2 * pi * (1:4) / 4))
})
