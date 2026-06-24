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
               unname(round(s$fit$par[res$parameter], 5L)))
})

test_that("bootstrap_rl and bootstrap_coef input validations", {
  s <- make_boot_fit()
  expect_error(bootstrap_rl(list(), data.frame()), regexp = "nhpp_fit")
  expect_error(bootstrap_coef(list(), data.frame()), regexp = "nhpp_fit")
  expect_error(bootstrap_rl(s$fit, list(), TRs = 10, R = 2), regexp = "must be a data frame")
  expect_error(bootstrap_coef(s$fit, list(), R = 2), regexp = "must be a data frame")
  expect_error(bootstrap_rl(s$fit, s$df, TRs = 10, R = 2, approach = "Z"), regexp = "must be one of 'A', 'B', 'C'")

  df_few <- s$df
  df_few$y <- 0
  expect_error(bootstrap_rl(s$fit, df_few, TRs = 10, R = 2), regexp = "fewer than 5 exceedances")
  expect_error(bootstrap_coef(s$fit, df_few, R = 2), regexp = "fewer than 5 exceedances")
})

test_that("bootstrap_rl and bootstrap_coef print verbose progress messages", {
  s <- make_boot_fit()
  expect_message(bootstrap_rl(s$fit, s$df, TRs = 10, R = 1L, approach = "A", verbose = TRUE), regexp = "replicate 1")
  expect_message(bootstrap_coef(s$fit, s$df, R = 1L, verbose = TRUE), regexp = "replicate 1")
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

  expect_error(
    bootstrap_coef(s$fit, s$df, R = 2, verbose = FALSE),
    regexp = "invalid GPD scale"
  )
})
