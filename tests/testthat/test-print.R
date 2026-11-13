# Tests for the print and plot methods of the result objects (print.R).

test_that("horizon_risk() returns a printable and plottable nhpp_horizon", {
  s <- make_bt_df()
  traj <- rep(list(data.frame(x = rep(0.5, s$n_obs))), 8L)
  hr <- horizon_risk(s$fit, c(12, 14), traj)
  expect_s3_class(hr, c("nhpp_horizon", "data.frame"), exact = TRUE)
  expect_identical(attr(hr, "n_trajectories"), 1L)
  out <- capture.output(print(hr))
  expect_match(out[1L], "horizon of 8 years \\(1 trajectory\\)")
  expect_false(any(grepl("cumulative", out)))
  capture.output(v <- withVisible(print(hr)))
  expect_false(v$visible)
  # column access and subsetting behave as for a data frame
  expect_type(hr$prob_exceed, "double")
  expect_s3_class(hr[1L, ], "data.frame")
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot(hr))
  expect_null(plot(hr[2L, ]))          # rows keep their curves
  bad <- hr
  attr(bad, "cumulative") <- NULL
  expect_error(plot(bad), "no cumulative")
})

test_that("print.nhpp_boot shows both tables and the bootstrap settings", {
  rl <- data.frame(TR = c(10, 50), RL_est = c(15.2, 18.9),
                   CI_low = c(14, 17), CI_high = c(16.5, 21), n_ok = 30L)
  cf <- data.frame(parameter = c("mu0", "mu.x"), estimate = c(5, 0.8),
                   CI_low = c(4.8, 0.6), CI_high = c(5.2, 1), p_zero = 0,
                   n_ok = 30L)
  b <- structure(list(rl = rl, coef = cf, rl_replicates = NULL,
                      coef_replicates = NULL),
                 info = list(R = 30L, level = 0.9, approach = "C",
                             occurrence = "fixed"),
                 class = "nhpp_boot")
  out <- capture.output(print(b))
  expect_match(out[1L], "R = 30 \\| 90% percentile intervals \\| occurrence: fixed")
  expect_true(any(grepl("approach C", out)))
  expect_true(any(grepl("mu.x", out)))
  b$rl <- NULL
  expect_false(any(grepl("Return levels", capture.output(print(b)))))
})

test_that("bootstrap_nhpp() returns an nhpp_boot object", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
                  obs_per_year = 50, verbose = FALSE)
  b <- bootstrap_nhpp(fit, df, TRs = NULL, R = 3L, verbose = FALSE)
  expect_s3_class(b, "nhpp_boot")
  expect_identical(attr(b, "info")$R, 3L)
  expect_true(is.data.frame(b$coef))
})

test_that("print.nhpp_backtest hides the legacy alias columns", {
  bt <- structure(list(
    results = data.frame(year = 2011:2015),
    binom_tests = data.frame(approach = "C", TR = 2, expected = 0.5,
                             observed = 0.4, n_exceedances = 2L, n_total = 5L,
                             p_value = 1, not_rejected = TRUE,
                             calibrated = TRUE),
    ks_test = NULL), class = "nhpp_backtest")
  out <- capture.output(print(bt))
  expect_match(out[1L], "5 validation years")
  expect_true(any(grepl("not_rejected", out)))
  expect_false(any(grepl("calibrated", out)))
  expect_false(any(grepl("Kolmogorov", out)))
  bt$ks_test <- data.frame(approach = "C", n = 5L, D_statistic = 0.2,
                           p_value = 0.9, not_rejected = TRUE,
                           calibrated = TRUE)
  out <- capture.output(print(bt))
  expect_true(any(grepl("Kolmogorov", out)))
  expect_false(any(grepl("calibrated", out)))
})

test_that("print.nhpp_var_generator reports deterministic and derived parts", {
  gen <- structure(list(p_opt = 1L, vars = c("x", "h"),
                        stochastic_vars = "x", n_harmonics = 1L,
                        period = 365, K = 1L, spectral_radius = 0.5,
                        derive = function(d) d, derive_lead = 3L),
                   class = "nhpp_var_generator")
  out <- capture.output(print(gen))
  expect_true(any(grepl("Deterministic   : h", out)))
  expect_true(any(grepl("1 harmonic,", out)))
  expect_true(any(grepl("lead of 3 observations", out)))
})

test_that("print.nhpp_var_generator summarises the generator", {
  set.seed(1)
  n <- 3650
  t <- rep(1:365, 10)
  df <- cbind(data.frame(x = as.numeric(stats::arima.sim(list(ar = 0.6), n))),
              seasonal_harmonics(t))
  df$y <- 10 + df$x + rexp(n, 1 / 3)
  fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
                  obs_per_year = 365, verbose = FALSE)
  gen <- fit_var_generator(fit, df)
  out <- capture.output(print(gen))
  expect_match(out[1L], sprintf("^VAR\\(%d\\) covariate generator", gen$p_opt))
  expect_true(any(grepl("Spectral radius", out)))
  expect_lt(length(out), 10L)
})
