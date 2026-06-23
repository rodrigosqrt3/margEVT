# Helper: enough data for walk-forward backtesting
make_backtest_df <- function(seed = 1L) {
  set.seed(seed)
  # 30 years of daily-ish data, 50 obs per year for speed
  n_years <- 30L
  n_per_year <- 50L
  n <- n_years * n_per_year
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20)),
    x    = stats::rnorm(n + 30L),
    year = rep(1991:2020, each = n_per_year)[seq_len(n + 30L)]
  )
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "none", verbose = FALSE)
  list(fit = fit, df = df)
}

test_that("backtest returns correct list structure", {
  s   <- make_backtest_df()
  res <- backtest(s$fit, s$df, varname = "y",
                  TRs = c(2, 5), n_obs = 50L,
                  window_years = 5L, min_train_years = 10L,
                  n_boot = 20L, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("results", "binom_tests", "ks_test"))
})

test_that("results data frame has correct columns", {
  s   <- make_backtest_df()
  res <- backtest(s$fit, s$df, varname = "y",
                  TRs = c(2, 5), n_obs = 50L,
                  window_years = 5L, min_train_years = 10L,
                  n_boot = 20L, verbose = FALSE)
  expect_true(all(c("window", "year", "M_obs", "F_ann",
                    "has_event", "exc_T2", "exc_T5") %in%
                    names(res$results)))
})

test_that("binom_tests has one row per TR", {
  s   <- make_backtest_df()
  res <- backtest(s$fit, s$df, varname = "y",
                  TRs = c(2, 5, 10), n_obs = 50L,
                  window_years = 5L, min_train_years = 10L,
                  n_boot = 20L, verbose = FALSE)
  expect_equal(nrow(res$binom_tests), 3L)
  expect_equal(res$binom_tests$TR, c(2, 5, 10))
})

test_that("F_ann values are in [0, 1]", {
  s   <- make_backtest_df()
  res <- backtest(s$fit, s$df, varname = "y",
                  TRs = 2, n_obs = 50L,
                  window_years = 5L, min_train_years = 10L,
                  n_boot = 20L, verbose = FALSE)
  f_vals <- res$results$F_ann[!is.na(res$results$F_ann)]
  expect_true(all(f_vals >= 0 & f_vals <= 1))
})

test_that("validation years are after training years in each window", {
  s   <- make_backtest_df()
  res <- backtest(s$fit, s$df, varname = "y",
                  TRs = 2, n_obs = 50L,
                  window_years = 5L, min_train_years = 10L,
                  n_boot = 20L, verbose = FALSE)
  # All validation years should be after min_train_years
  expect_true(all(res$results$year >= 1991 + 10L))
})

test_that("missing varname throws informative error", {
  s <- make_backtest_df()
  expect_error(
    backtest(s$fit, s$df, varname = "nonexistent",
             n_obs = 50L, verbose = FALSE),
    regexp = "nonexistent"
  )
})

test_that("not enough years throws informative error", {
  s      <- make_backtest_df()
  df_sub <- s$df[s$df$year <= 1995, ]
  expect_error(
    backtest(s$fit, df_sub, varname = "y",
             n_obs = 50L, min_train_years = 20L, verbose = FALSE),
    regexp = "min_train_years"
  )
})

test_that("backtest rejects non nhpp_fit", {
  expect_error(backtest(list(), data.frame(), varname = "y"),
               regexp = "nhpp_fit")
})
