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

test_that("bootstrap_rl rejects non nhpp_fit", {
  expect_error(bootstrap_rl(list(), data.frame()),
               regexp = "nhpp_fit")
})

test_that("bootstrap_coef rejects non nhpp_fit", {
  expect_error(bootstrap_coef(list(), data.frame()),
               regexp = "nhpp_fit")
})
