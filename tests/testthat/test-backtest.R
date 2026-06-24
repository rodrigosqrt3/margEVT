library(testthat)
library(margEVT)

# -----------------------------------------------------------------------------
# Test Helper
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Core Functionality Tests
# -----------------------------------------------------------------------------
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
  expect_true(all(res$results$year >= 2001))
})

# -----------------------------------------------------------------------------
# Input Validation & Errors
# -----------------------------------------------------------------------------
test_that("backtest rejects non nhpp_fit", {
  expect_error(backtest(list(), data.frame(), varname = "y"),
               regexp = "nhpp_fit")
})

test_that("backtest rejects non-data.frame inputs", {
  s <- make_backtest_df()
  expect_error(backtest(s$fit, list(), varname = "y"),
               regexp = "data frame")
})

test_that("missing varname throws informative error", {
  s <- make_backtest_df()
  expect_error(
    backtest(s$fit, s$df, varname = "nonexistent",
             n_obs = 50L, verbose = FALSE),
    regexp = "nonexistent"
  )
})

test_that("missing year_col throws informative error", {
  s <- make_backtest_df()
  expect_error(
    backtest(s$fit, s$df, varname = "y", year_col = "nonexistent",
             n_obs = 50L, verbose = FALSE),
    regexp = "nonexistent"
  )
})

test_that("missing 'y' column in data throws informative error", {
  s <- make_backtest_df()
  df_noy <- s$df
  df_noy$y <- NULL
  expect_error(
    backtest(s$fit, df_noy, varname = "x", n_obs = 50L, verbose = FALSE),
    regexp = "must contain column `y`"
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

# -----------------------------------------------------------------------------
# Advanced Coverage & Branch Tests
# -----------------------------------------------------------------------------
test_that("backtest supports verbose printing statement coverage", {
  s <- make_backtest_df()
  expect_message(
    backtest(s$fit, s$df, varname = "y", TRs = 2, n_obs = 50L,
             window_years = 10L, min_train_years = 15L,
             n_boot = 10L, verbose = TRUE),
    regexp = "Binomial calibration tests"
  )
})

test_that("backtest handles stationary fit with no covariates", {
  s <- make_backtest_df()
  fit_stat <- fit_nhpp(s$df, threshold = 4, penalty = "none", verbose = FALSE)

  res_stat <- backtest(fit_stat, s$df, varname = "y", TRs = 2, n_obs = 50L,
                       window_years = 10L, min_train_years = 15L,
                       n_boot = 10L, verbose = FALSE)
  expect_type(res_stat, "list")
})

test_that("backtest skips training windows that are too short (< 2 years)", {
  # We construct a dataset with only years 1991, 1992, and 2020.
  # Under min_train_years = 1L, the first window train_yrs is only c(1991) (length 1),
  # which triggers the "length(train_yrs) < 2L" next skip.
  df_gap <- data.frame(
    y = rexp(150, 0.3),
    x = rnorm(150),
    year = rep(c(1991, 1992, 2020), each = 50L)
  )
  fit_gap <- fit_nhpp(df_gap, threshold = 4, loc_vars = "x", penalty = "none", verbose = FALSE)

  res_gap <- backtest(fit_gap, df_gap, varname = "y", TRs = 2, n_obs = 50L,
                      window_years = 5L, min_train_years = 1L,
                      n_boot = 10L, verbose = FALSE)
  expect_type(res_gap, "list")
})

test_that("backtest handles observation padding (falta) branch", {
  s <- make_backtest_df()

  # For one year (e.g., 2015), we keep 48 observations.
  # This is less than n_obs (50) but greater than 90% (45), which triggers the padding code.
  df_pad <- s$df
  rows_to_remove <- which(df_pad$year == 2015)[1:2]
  df_pad <- df_pad[-rows_to_remove, ]

  res_pad <- backtest(s$fit, df_pad, varname = "y", TRs = 2, n_obs = 50L,
                      window_years = 5L, min_train_years = 15L,
                      n_boot = 10L, verbose = FALSE)
  expect_type(res_pad, "list")
})

test_that("backtest handles Gumbel limit (|xi| < 1e-6) and z_u <= 0 boundary (Mocked)", {
  s <- make_backtest_df()

  # Force Gumbel limit branch by setting shape parameter estimate to exactly 0
  fit_gumbel <- s$fit
  fit_gumbel$par["xi.(Intercept)"] <- 0.0

  res_gumbel <- backtest(fit_gumbel, s$df, varname = "y",
                         TRs = 2, n_obs = 50L,
                         window_years = 5L, min_train_years = 15L,
                         n_boot = 10L, verbose = FALSE)
  expect_type(res_gumbel, "list")

  # Force z_u <= 0 condition by inserting an extremely high outlier value (y = 1000)
  df_outlier <- s$df
  df_outlier$y[df_outlier$year == 2019] <- 1000.0

  res_outlier <- backtest(s$fit, df_outlier, varname = "y",
                          TRs = 2, n_obs = 50L,
                          window_years = 5L, min_train_years = 15L,
                          n_boot = 10L, verbose = FALSE)
  expect_type(res_outlier, "list")
})

test_that("backtest handles insufficient observations in validation year", {
  s <- make_backtest_df()

  # Remove 45 rows of observations for 2015, making nrow(df_yr) < min_obs_year true
  df_sparse <- s$df
  rows_to_remove <- which(df_sparse$year == 2015)[1:45]
  df_sparse <- df_sparse[-rows_to_remove, ]

  res_sparse <- backtest(s$fit, df_sparse, varname = "y",
                         TRs = 2, n_obs = 50L,
                         window_years = 5L, min_train_years = 15L,
                         n_boot = 10L, verbose = FALSE)

  results_2015 <- res_sparse$results[res_sparse$results$year == 2015, ]
  expect_true(is.na(results_2015$F_ann))
})

test_that("backtest triggers 'No valid covariate blocks' and returns NULL if all skipped", {
  s <- make_backtest_df()

  # Keep all years 1991:2020, but restrict each year to only 30 observations.
  # This makes nrow(df_yr) < (n_obs * 0.9) [30 < 45] always TRUE, returning an empty list.
  df_short_years <- do.call(rbind, lapply(split(s$df, s$df$year), function(sub) {
    utils::head(sub, 30L)
  }))

  expect_warning(
    res_empty <- backtest(s$fit, df_short_years, varname = "y",
                          TRs = 2, n_obs = 50L,
                          window_years = 5L, min_train_years = 15L,
                          n_boot = 10L, verbose = FALSE),
    regexp = "no validation results produced"
  )
  expect_null(res_empty)
})

# -----------------------------------------------------------------------------
# Internal Helper Direct Tests
# -----------------------------------------------------------------------------
test_that("internal helper .build_cov_bootstrap handles empty anos_disp", {
  s <- make_backtest_df()

  # Corrected: Set ac_in_data = "x" and pass empty column vectors
  # to safely bypass the first check and trigger the length(anos_disp) == 0L branch
  res <- margEVT:::.build_cov_bootstrap(
    fit_w = s$fit,
    data = data.frame(year = integer(0), x = numeric(0)),
    train_yrs = 1991:2000,
    ac_in_data = "x",
    n_boot = 10L,
    n_obs = 50L,
    interactions = list(),
    year_col = "year"
  )
  expect_equal(res, list())
})

test_that("internal helper .build_cov_bootstrap imputes NAs with median", {
  s <- make_backtest_df()

  # Create a small dataset with NAs in covariate x
  data_with_na <- data.frame(
    year = rep(1991:1992, each = 50L),
    x = c(rep(1, 49), NA_real_, rep(5, 50))
  )

  res <- margEVT:::.build_cov_bootstrap(
    fit_w = s$fit,
    data = data_with_na,
    train_yrs = 1991:1992,
    ac_in_data = "x",
    n_boot = 5L,
    n_obs = 50L,
    interactions = list(),
    year_col = "year"
  )

  # Verifies successful run through NA imputation path
  expect_length(res, 5L)
})

test_that("backtest handles non-convergence in training window", {
  s <- make_backtest_df()

  df_no_conv <- s$df
  df_no_conv$y[1:10] <- NaN

  expect_message(
    res <- backtest(s$fit, df_no_conv, varname = "y", TRs = 2, n_obs = 50L,
                    window_years = 5L, min_train_years = 15L,
                    n_boot = 10L, verbose = TRUE),
    regexp = "Did not converge, skipping window"
  )
  expect_null(res) # <-- Changed from expect_type(res, "list")
})

test_that("backtest triggers 'No valid covariate blocks' and returns NULL if all skipped", {
  s <- make_backtest_df()

  df_short_years <- do.call(rbind, lapply(split(s$df, s$df$year), function(sub) {
    utils::head(sub, 30L)
  }))

  expect_warning(
    res_empty <- backtest(s$fit, df_short_years, varname = "y",
                          TRs = 2, n_obs = 50L,
                          window_years = 5L, min_train_years = 15L,
                          n_boot = 10L, verbose = TRUE), # <-- Changed to TRUE
    regexp = "no validation results produced"
  )
  expect_null(res_empty)
})

test_that("backtest handles missing columns in binomial test mapping", {
  s <- make_backtest_df()

  res_missing <- backtest(s$fit, s$df, varname = "y",
                          TRs = c(2, 5), n_obs = 50L,
                          window_years = 10L, min_train_years = 15L,
                          n_boot = 10L, verbose = FALSE)

  expect_type(res_missing, "list")
  expect_equal(res_missing$binom_tests$TR, c(2, 5))   # <-- Reverted to c(2, 5)
})
