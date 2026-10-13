# Tests for the 0.4.0 bootstrap and backtest machinery: warm-started refits
# with the original penalty weights, per-replicate random-number streams,
# parallel execution, trend renormalization, and fixed covariates.

.boot_data <- function(seed = 1L, n_years = 30L, n_per_year = 50L) {
  set.seed(seed)
  n <- n_years * n_per_year
  df <- data.frame(
    x    = stats::rnorm(n),
    z    = stats::rnorm(n),
    year = rep(seq_len(n_years) + 1990L, each = n_per_year)
  )
  df$trend <- (df$year - min(df$year)) / (max(df$year) - min(df$year))
  df$y <- 5 + df$x + 2 * df$trend + stats::rexp(n, 1 / 2)
  df
}

# TRUE when the installed margEVT is this version, so PSOCK workers (which
# load the installed package) run the code under test.
.installed_matches <- function() {
  ip <- tryCatch(utils::installed.packages()[, "Version"],
                 error = function(e) character(0L))
  "margEVT" %in% names(ip) &&
    identical(unname(ip["margEVT"]),
              as.character(utils::packageVersion("margEVT")))
}

test_that("RNG streams are reproducible, independent and seed-preserving", {
  set.seed(3L); before <- stats::runif(2L)
  set.seed(3L); first <- stats::runif(1L)
  s1 <- margEVT:::.rng_streams(3L, 10L)
  second <- stats::runif(1L)
  expect_equal(c(first, second), before)
  expect_identical(s1, margEVT:::.rng_streams(3L, 10L))
  expect_false(identical(s1[[1L]], s1[[2L]]))

  u1 <- margEVT:::.with_rng_stream(s1[[2L]], stats::runif(3L))
  u2 <- margEVT:::.with_rng_stream(s1[[2L]], stats::runif(3L))
  expect_identical(u1, u2)
  expect_identical(RNGkind()[1L], "Mersenne-Twister")

  # seed = NULL draws the base seed from (and advances) the caller stream
  set.seed(4L); a <- margEVT:::.rng_streams(1L)
  set.seed(4L); b <- margEVT:::.rng_streams(1L)
  expect_identical(a, b)

  # without a pre-existing seed none is left behind
  old <- get(".Random.seed", envir = .GlobalEnv)
  rm(".Random.seed", envir = .GlobalEnv)
  margEVT:::.rng_streams(1L, 1L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", old, envir = .GlobalEnv)
})

test_that("bootstrap is reproducible, keeps the caller RNG and stores replicates", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
                  obs_per_year = 50, verbose = FALSE)
  set.seed(8L); expected <- stats::runif(2L)
  set.seed(8L); first <- stats::runif(1L)
  a <- bootstrap_coef(fit, df, R = 6L, seed = 5L, verbose = FALSE)
  second <- stats::runif(1L)
  expect_equal(c(first, second), expected)
  b <- bootstrap_coef(fit, df, R = 6L, seed = 5L, verbose = FALSE)
  expect_identical(a, b)
  expect_equal(dim(attr(a, "replicates")), c(6L, length(fit$par)))
  expect_identical(colnames(attr(a, "replicates")), names(fit$par))

  rl <- bootstrap_rl(fit, df, TRs = c(5, 10), R = 3L, approach = "A",
                     seed = NULL, verbose = FALSE)
  expect_identical(colnames(attr(rl, "replicates")), c("T5", "T10"))
})

test_that("bootstrap validates seed, cores and approach", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
                  obs_per_year = 50, verbose = FALSE)
  expect_error(bootstrap_coef(fit, df, R = 2L, seed = "a"), "seed")
  expect_error(bootstrap_coef(fit, df, R = 2L, cores = 0), "cores")
  expect_error(bootstrap_rl(fit, df, TRs = 10, R = 2L, cores = 0), "cores")
  expect_error(bootstrap_rl(fit, df, TRs = 10, R = 2L, approach = c("A", "B")),
               "approach")
})

test_that("bootstrap refits keep the original penalty weights and settings", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = c("x", "z"),
                  penalty = "lasso", lambda = 2, penalty_weights = "adaptive",
                  grad_tol = 1e-3, adaptive_gamma = 0.5,
                  obs_per_year = 50, verbose = FALSE)
  expect_identical(fit$grad_tol, 1e-3)
  expect_identical(fit$adaptive_gamma, 0.5)
  fit_b <- margEVT:::.refit_boot(df, fit)
  expect_identical(fit_b$penalty_weights, fit$penalty_weights)
  expect_identical(fit_b$weights_source, "user")
  expect_identical(fit_b$grad_tol, 1e-3)
  expect_equal(fit_b$par, fit$par, tolerance = 1e-4)

  fit_a <- margEVT:::.refit_boot(df[df$year < 2010, ], fit,
                                 penalty_weights = "adaptive")
  expect_match(fit_a$weights_source, "^adaptive")
  expect_identical(fit_a$adaptive_gamma, 0.5)

  # legacy objects without the 0.4.0 fields still refit
  legacy <- fit
  legacy[c("penalty_weights", "lambda_scaling", "adaptive_gamma",
           "grad_tol", "active_tol")] <- list(NULL)
  expect_s3_class(margEVT:::.refit_boot(df, legacy), "nhpp_fit")
})

test_that("parallel bootstrap and backtest equal their serial versions", {
  skip_on_cran()
  skip_if_not(.installed_matches(),
              "installed margEVT differs from the tested version")
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
                  obs_per_year = 50, verbose = FALSE)
  expect_identical(
    bootstrap_coef(fit, df, R = 4L, seed = 2L, cores = 2L, verbose = FALSE),
    bootstrap_coef(fit, df, R = 4L, seed = 2L, cores = 1L, verbose = FALSE))
  args <- list(fit, df, varname = "y", TRs = 2, n_obs = 50L,
               window_years = 10L, min_train_years = 15L, n_boot = 10L,
               verbose = FALSE)
  expect_identical(do.call(backtest, c(args, cores = 2L)),
                   do.call(backtest, c(args, cores = 1L)))
})

test_that("backtest renormalizes a trend covariate within each window", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = c("x", "trend"),
                  penalty = "none", obs_per_year = 50, verbose = FALSE)
  bt <- backtest(fit, df, varname = "y", TRs = c(2, 5), n_obs = 50L,
                 window_years = 5L, min_train_years = 15L, n_boot = 20L,
                 trend_col = "trend", verbose = FALSE)
  bt0 <- backtest(fit, df, varname = "y", TRs = c(2, 5), n_obs = 50L,
                  window_years = 5L, min_train_years = 15L, n_boot = 20L,
                  verbose = FALSE)
  expect_true(all(bt$results$F_ann >= 0 & bt$results$F_ann <= 1,
                  na.rm = TRUE))
  # the window-specific trend changes the validation PIT values
  expect_false(isTRUE(all.equal(bt$results$F_ann, bt0$results$F_ann)))

  # the trend is also propagated to interactions that involve it
  df$xt <- df$x * df$trend
  fit_i <- fit_nhpp(df, threshold = 10, loc_vars = c("x", "trend", "xt"),
                    penalty = "none", obs_per_year = 50, verbose = FALSE)
  bt_i <- backtest(fit_i, df, varname = "y", TRs = 2, n_obs = 50L,
                   window_years = 10L, min_train_years = 15L, n_boot = 10L,
                   trend_col = "trend",
                   interactions = list(xt = c("x", "trend")),
                   verbose = FALSE)
  expect_s3_class(bt_i$results, "data.frame")

  expect_error(backtest(fit, df, "y", trend_col = "nope"), "trend_col")
  expect_error(backtest(fit, df, "y", trend_col = c("x", "trend")),
               "trend_col")
})

test_that("backtest accepts fixed covariates, BIC re-selection and adaptive fits", {
  df <- .boot_data()
  fit <- fit_nhpp(df, threshold = 10, loc_vars = c("x", "z"),
                  penalty = "lasso", lambda = 1, penalty_weights = "adaptive",
                  obs_per_year = 50, verbose = FALSE)
  bt <- backtest(fit, df, varname = "y", TRs = 2, n_obs = 50L,
                 window_years = 10L, min_train_years = 15L, n_boot = 10L,
                 fixed = list(z = 0), reselect_lambda = TRUE,
                 verbose = FALSE)
  expect_s3_class(bt$binom_tests, "data.frame")
  # covariates inactive in a window refit are zero-filled without warnings
  fit_z <- fit_nhpp(df, threshold = 10, loc_vars = c("x", "z"),
                    penalty = "lasso", lambda = 1, penalty_weights = "adaptive",
                    obs_per_year = 50, verbose = FALSE)
  expect_false("z" %in% active_covariates(fit_z))
  expect_warning(
    backtest(fit_z, df, varname = "y", TRs = 2, n_obs = 50L,
             window_years = 10L, min_train_years = 15L, n_boot = 10L,
             verbose = FALSE),
    NA)
  expect_error(backtest(fit, df, "y", fixed = list(1)), "fixed")
  expect_error(backtest(fit, df, "y", reselect_lambda = NA), "reselect_lambda")
  expect_error(backtest(fit, df, "y", cores = 0), "cores")
})
