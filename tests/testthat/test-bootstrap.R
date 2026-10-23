library(testthat)
library(margEVT)

# Helper: fit a simple model we can bootstrap from
make_boot_fit <- function(seed = 1L) {
  set.seed(seed)
  n  <- 400L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(40L, 5, 20)),
    x    = stats::rnorm(n + 40L),
    year = rep(1981:2021, each = 365L)[seq_len(n + 40L)]
  )
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "none", verbose = FALSE)
  list(fit = fit, df = df)
}

test_that("bootstrap_rl returns correct structure", {
  s   <- make_boot_fit()
  res <- bootstrap_rl(s$fit, s$df,
                      TRs      = c(10, 50),
                      R        = 10L,
                      approach = "A",
                      verbose  = FALSE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("TR", "RL_est", "CI_low", "CI_high", "n_ok"))
  expect_equal(nrow(res), 2L)
})

test_that("bootstrap_rl CI_low < RL_est < CI_high", {
  s   <- make_boot_fit()
  res <- bootstrap_rl(s$fit, s$df,
                      TRs      = c(10, 50),
                      R        = 20L,
                      approach = "A",
                      verbose  = FALSE)
  expect_true(all(res$CI_low  < res$RL_est,  na.rm = TRUE))
  expect_true(all(res$RL_est  < res$CI_high, na.rm = TRUE))
})

test_that("bootstrap_rl n_ok <= R", {
  s   <- make_boot_fit()
  res <- bootstrap_rl(s$fit, s$df,
                      TRs      = 10,
                      R        = 15L,
                      approach = "A",
                      verbose  = FALSE)
  expect_lte(res$n_ok, 15L)
})

test_that("bootstrap_rl approach C works with year column", {
  s   <- make_boot_fit()
  res <- bootstrap_rl(s$fit, s$df,
                      TRs              = 10,
                      R                = 10L,
                      approach         = "C",
                      marginalize_args = list(n_boot = 20L, n_obs = 11L, seed = 1L),
                      verbose          = FALSE)
  expect_true(is.finite(res$RL_est))
})

test_that("bootstrap_coef returns correct structure", {
  s   <- make_boot_fit()
  res <- bootstrap_coef(s$fit, s$df, R = 10L, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("parameter", "estimate", "CI_low", "CI_high", "n_ok"))
  expect_equal(nrow(res), length(s$fit$par))
})

test_that("bootstrap_coef estimates match fit$par", {
  s   <- make_boot_fit()
  res <- bootstrap_coef(s$fit, s$df, R = 10L, verbose = FALSE)
  expect_equal(res$estimate,
               unname(round(s$fit$par_oper[res$parameter], 5L)))
})

test_that("bootstrap refits preserve the original design and free variables", {
  set.seed(11L)
  n <- 500L
  tt <- seq_len(n)
  df <- data.frame(
    y = c(stats::rexp(n - 40L, 0.3), stats::runif(40L, 5, 20)),
    x = stats::rnorm(n),
    cos1 = cos(2 * pi * tt / 100),
    sen1 = sin(2 * pi * tt / 100)
  )
  fit <- fit_nhpp(
    df, threshold = 4,
    loc_vars = c("cos1", "sen1", "x"),
    scale_vars = c("cos1", "x"),
    free_vars = c("cos1", "sen1"),
    penalty = "lasso", lambda = 0.5,
    obs_per_year = 100, verbose = FALSE
  )

  fit_b <- margEVT:::.refit_boot(df, fit)
  expect_s3_class(fit_b, "nhpp_fit")
  expect_identical(colnames(fit_b$dm$X_mu), colnames(fit$dm$X_mu))
  expect_identical(colnames(fit_b$dm$X_sigma), colnames(fit$dm$X_sigma))
  expect_identical(fit_b$dm$idx_pen_mu, fit$dm$idx_pen_mu)
  expect_identical(fit_b$dm$idx_pen_sigma, fit$dm$idx_pen_sigma)
})

test_that("bootstrap_rl and bootstrap_coef input validations", {
  s <- make_boot_fit()
  expect_error(bootstrap_rl(list(), data.frame()), regexp = "nhpp_fit")
  expect_error(bootstrap_coef(list(), data.frame()), regexp = "nhpp_fit")
  expect_error(bootstrap_rl(s$fit, list(), TRs = 10, R = 2), regexp = "must be a data frame")
  expect_error(bootstrap_coef(s$fit, list(), R = 2), regexp = "must be a data frame")
  expect_error(bootstrap_rl(s$fit, s$df, TRs = 10, R = 2, approach = "Z"), regexp = "must be one of 'A', 'B', 'C'")
  expect_error(bootstrap_rl(s$fit, s$df, TRs = 1, R = 2), regexp = "greater than 1")
  expect_error(bootstrap_rl(s$fit, s$df, TRs = 10, R = 0), regexp = "positive integer")
  expect_error(bootstrap_rl(s$fit, s$df, TRs = 10, R = 2, level = 1), regexp = "between 0 and 1")
  expect_error(bootstrap_coef(s$fit, s$df, R = 0), regexp = "positive integer")
  expect_error(bootstrap_coef(s$fit, s$df, R = 2, level = 1), regexp = "between 0 and 1")

  df_few <- s$df
  df_few$y <- 0
  expect_error(bootstrap_rl(s$fit, df_few, TRs = 10, R = 2), regexp = "fewer than 5 exceedances")
  expect_error(bootstrap_coef(s$fit, df_few, R = 2), regexp = "fewer than 5 exceedances")
})

test_that("bootstrap functions reject non-finite response values", {
  s <- make_boot_fit()
  bad <- s$df
  bad$y[1L] <- NA_real_

  expect_error(
    bootstrap_rl(s$fit, bad, TRs = 10, R = 2L, verbose = FALSE),
    regexp = "data\\$y.*finite numeric values"
  )
  expect_error(
    bootstrap_coef(s$fit, bad, R = 2L, verbose = FALSE),
    regexp = "data\\$y.*finite numeric values"
  )
})

test_that("GPD simulation rejects non-finite generated excesses", {
  expect_error(
    margEVT:::.simulate_exceedances(
      mu_t = 4, sigma_t = 1, xi_t = 1e308,
      threshold = 4, exc_idx = 1L
    ),
    regexp = "generated an invalid GPD excess"
  )
})

test_that("bootstrap_rl and bootstrap_coef print verbose progress messages", {
  s <- make_boot_fit()
  expect_message(bootstrap_rl(s$fit, s$df, TRs = 10, R = 1L, approach = "A", verbose = TRUE), regexp = "replicates succeeded")
  expect_message(bootstrap_coef(s$fit, s$df, R = 1L, verbose = TRUE), regexp = "replicates succeeded")
})

test_that("bootstrap functions handle non-convergence skips", {
  s <- make_boot_fit()
  mock_fit_nhpp <- function(...) list(converged = FALSE)
  testthat::local_mocked_bindings(fit_nhpp = mock_fit_nhpp)

  res_rl <- bootstrap_rl(s$fit, s$df, TRs = 10, R = 1L, approach = "A", verbose = FALSE)
  expect_equal(res_rl$n_ok, 0L)
  res_coef <- bootstrap_coef(s$fit, s$df, R = 1L, verbose = FALSE)
  expect_equal(res_coef$n_ok[1], 0L)
})

test_that("bootstrap_rl handles NULL from marginalize", {
  s <- make_boot_fit()
  mock_marg <- function(...) NULL
  testthat::local_mocked_bindings(marginalize = mock_marg)

  res <- bootstrap_rl(s$fit, s$df, TRs = 10, R = 1L, approach = "A", verbose = FALSE)
  expect_equal(res$n_ok, 0L)
})

test_that("bootstrap catches invalid GPD scale limit", {
  s <- make_boot_fit()

  # Directly overwrite the stored fitted values!
  # beta_t = sigma + xi * (threshold - mu)
  # beta_t = exp(-50) + (-10) * (4 - 0) = -40 <= 0
  s$fit$fitted$mu    <- rep(0.0, nrow(s$df))
  s$fit$fitted$sigma <- rep(exp(-50.0), nrow(s$df))
  s$fit$fitted$xi    <- rep(-10.0, nrow(s$df))
  s$fit$fitted_oper <- s$fit$fitted

  expect_error(
    bootstrap_coef(s$fit, s$df, R = 2, verbose = FALSE),
    regexp = "invalid GPD scale"
  )
})

# -----------------------------------------------------------------------------
# Reproducibility, refit settings and parallel execution
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# refit_nhpp(), bootstrap_nhpp() and simulated occurrences
# -----------------------------------------------------------------------------

test_that("refit_nhpp refits with the same specification", {
  s <- make_bt_df()
  d2 <- s$df; set.seed(9); d2$y <- 5 + 0.8 * d2$x + stats::rexp(nrow(d2), 1 / 2)
  f2 <- refit_nhpp(s$fit, d2)
  expect_s3_class(f2, "nhpp_fit")
  expect_true(f2$converged)
  expect_identical(names(f2$par), names(s$fit$par))
  expect_error(refit_nhpp(list(), d2), "nhpp_fit")
  expect_error(refit_nhpp(s$fit, d2[, "x", drop = FALSE]), "column `y`")
})

test_that("bootstrap_nhpp reproduces bootstrap_rl and bootstrap_coef", {
  s <- make_bt_df()
  b  <- bootstrap_nhpp(s$fit, s$df, TRs = c(5, 20), R = 12, approach = "C",
                       marginalize_args = list(year_col = "year", n_boot = 20),
                       seed = 7, verbose = FALSE)
  br <- bootstrap_rl(s$fit, s$df, TRs = c(5, 20), R = 12, approach = "C",
                     marginalize_args = list(year_col = "year", n_boot = 20),
                     seed = 7, verbose = FALSE)
  bc <- bootstrap_coef(s$fit, s$df, R = 12, seed = 7, verbose = FALSE)
  expect_equal(b$rl$CI_low, br$CI_low)
  expect_equal(b$rl$CI_high, br$CI_high)
  expect_equal(b$coef$CI_low, bc$CI_low)
  expect_equal(dim(b$coef_replicates), c(12L, length(s$fit$par)))
  expect_null(bootstrap_nhpp(s$fit, s$df, TRs = NULL, R = 5, seed = 1,
                             verbose = FALSE)$rl)
})

test_that("bootstrap_nhpp validates inputs, reports progress and keeps failures", {
  s <- make_bt_df()
  expect_error(bootstrap_nhpp(s$fit, s$df, R = 2, marginalize_args = 1,
                              verbose = FALSE), "marginalize_args")
  expect_error(bootstrap_nhpp(s$fit, s$df, R = 2, approach = "D",
                              verbose = FALSE), "approach")
  expect_error(bootstrap_nhpp(s$fit, s$df, TRs = 1, R = 2, verbose = FALSE),
               "TRs")
  # approach B without `mc_sample`: marginalize() fails, levels become NA
  msgs <- testthat::capture_messages(
    b <- bootstrap_nhpp(s$fit, s$df, TRs = 10, R = 2, approach = "B",
                        seed = 1, verbose = TRUE)
  )
  expect_match(msgs, "bootstrap_nhpp: R=2 \\| approach=B", all = FALSE)
  expect_match(msgs, "replicates succeeded", all = FALSE)
  expect_true(is.na(b$rl$RL_est))
  expect_true(all(is.na(b$rl_replicates)))
})

test_that("bootstrap can resample the occurrence pattern", {
  s <- make_bt_df()
  bf <- bootstrap_coef(s$fit, s$df, R = 8, seed = 3, verbose = FALSE)
  bs <- bootstrap_coef(s$fit, s$df, R = 8, seed = 3, occurrence = "simulated",
                       verbose = FALSE)
  expect_true(all(bs$n_ok > 0))
  rf <- attr(bf, "replicates"); rs <- attr(bs, "replicates")
  expect_false(isTRUE(all.equal(rf, rs)))           # different resampling scheme
  expect_identical(bootstrap_coef(s$fit, s$df, R = 8, seed = 3,
                                  occurrence = "fixed", verbose = FALSE), bf)
  expect_error(bootstrap_coef(s$fit, s$df, R = 2, occurrence = "x",
                              verbose = FALSE), "should be one of")
})

test_that("simulated occurrences require the threshold inside the support", {
  s <- make_bt_df()
  testthat::local_mocked_bindings(.tail_measure_at_level = function(...) Inf)
  expect_error(bootstrap_coef(s$fit, s$df, R = 2, occurrence = "simulated",
                              verbose = FALSE),
               "threshold lies outside the fitted support")
})

test_that("replicates with fewer than 5 simulated exceedances are dropped", {
  s <- make_bt_df()
  testthat::local_mocked_bindings(.tail_measure_at_level = function(...) 0)
  b <- bootstrap_nhpp(s$fit, s$df, TRs = NULL, R = 3, seed = 1,
                      occurrence = "simulated", verbose = FALSE)
  expect_true(all(is.na(b$coef_replicates)))
  expect_true(all(b$coef$n_ok == 0L))
})
