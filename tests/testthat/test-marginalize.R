library(testthat)
library(margEVT)

make_stationary_fit <- function(seed = 1L) {
  set.seed(seed)
  n  <- 500L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(40, 5, 20)),
    year = rep(1981:2020, each = 13L)[seq_len(n + 40L)]
  )
  list(
    fit = fit_nhpp(df, threshold = 4, lambda = 0, verbose = FALSE),
    df  = df
  )
}

make_cov_fit <- function(seed = 2L) {
  set.seed(seed)
  n  <- 500L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(40, 5, 20)),
    x    = stats::rnorm(n + 40L),
    year = rep(1981:2020, each = 13L)[seq_len(n + 40L)]
  )
  list(
    fit = fit_nhpp(df, threshold = 4, loc_vars = "x",
                   lambda = 0, verbose = FALSE),
    df  = df
  )
}

test_that("approach A returns correct structure", {
  s   <- make_stationary_fit()
  res <- marginalize(s$fit, s$df, TRs = c(10, 50), approaches = "A")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("approach", "scenario", "TR", "RL") %in% names(res)))
  expect_equal(nrow(res), 2L)
  expect_true(all(res$approach == "A"))
})

test_that("approach A with multiple scenarios returns all rows", {
  s   <- make_stationary_fit()
  sc  <- list(low = list(), mid = list(), high = list())
  res <- marginalize(s$fit, s$df, TRs = c(10, 50),
                     approaches = "A", scenarios = sc)
  expect_equal(nrow(res), 6L)   # 3 scenarios x 2 TRs
  expect_equal(sort(unique(res$scenario)), c("high", "low", "mid"))
})

test_that("return levels increase with return period", {
  s   <- make_stationary_fit()
  res <- marginalize(s$fit, s$df,
                     TRs = c(2, 10, 50, 100), approaches = "A")
  expect_true(all(diff(res$RL) > 0))
})

test_that("tail measure respects lower and upper GEV endpoints", {
  lower_endpoint <- margEVT:::.tail_measure_at_level(
    z = -6, mu_t = 0, sigma_t = 1, xi_t = 0.2
  )
  upper_endpoint <- margEVT:::.tail_measure_at_level(
    z = 6, mu_t = 0, sigma_t = 1, xi_t = -0.2
  )

  expect_identical(lower_endpoint, Inf)
  expect_identical(upper_endpoint, 0)
})

test_that("annual probabilities use the correct endpoint convention", {
  below_lower <- margEVT:::.annual_exceedance_prob(
    z = -6, mu_t = 0, sigma_t = 1, xi_t = 0.2, n_obs = 1L
  )
  above_upper <- margEVT:::.annual_exceedance_prob(
    z = 6, mu_t = 0, sigma_t = 1, xi_t = -0.2, n_obs = 1L
  )

  expect_identical(below_lower, 0)
  expect_identical(above_upper, 1)
})

test_that("annual probabilities reject malformed trajectory stacks", {
  expect_error(
    margEVT:::.annual_exceedance_prob(
      z = 1, mu_t = 1:3, sigma_t = rep(1, 3), xi_t = rep(0, 3),
      n_obs = 2L
    ),
    regexp = "divisible"
  )
})

test_that("annual probabilities reject invalid fitted parameters", {
  expect_error(
    margEVT:::.annual_exceedance_prob(
      z = 1, mu_t = c(0, NA), sigma_t = c(1, 1),
      xi_t = c(0, 0), n_obs = 2L
    ),
    regexp = "finite"
  )
  expect_error(
    margEVT:::.annual_exceedance_prob(
      z = 1, mu_t = c(0, 0), sigma_t = c(1, 0),
      xi_t = c(0, 0), n_obs = 2L
    ),
    regexp = "strictly positive"
  )
})

test_that("return-level search does not extrapolate below the threshold", {
  root <- margEVT:::.find_return_level(
    TR = 2,
    f_annual = function(z) 0.9,
    z_lo = 10,
    z_hi = 20
  )
  expect_true(is.na(root))
})

test_that("return-level search expands only the upper bracket", {
  root <- margEVT:::.find_return_level(
    TR = 10,
    f_annual = stats::plogis,
    z_lo = 0,
    z_hi = 1
  )
  expect_equal(root, stats::qlogis(0.9), tolerance = 1e-4)
})

test_that("approach C returns correct structure", {
  s   <- make_stationary_fit()
  res <- marginalize(s$fit, s$df, TRs = c(10, 50),
                     approaches = "C", n_boot = 50L, seed = 1L)
  expect_equal(nrow(res), 2L)
  expect_true(all(res$approach == "C"))
  expect_true(all(is.finite(res$RL)))
})

test_that("approaches A and C give similar results for stationary model", {
  s    <- make_stationary_fit()
  resA <- marginalize(s$fit, s$df, TRs = 10, approaches = "A")
  resC <- marginalize(s$fit, s$df, TRs = 10, approaches = "C",
                      n_boot = 200L, seed = 1L)
  expect_lt(abs(resA$RL - resC$RL) / resA$RL, 0.20)
})

test_that("approach B with mc_sample returns correct structure", {
  s  <- make_cov_fit()
  set.seed(1L)
  mc <- lapply(seq_len(50L), function(i)
    data.frame(x = stats::rnorm(365L)))
  res <- marginalize(s$fit, s$df, TRs = c(10, 50),
                     approaches = "B", mc_sample = mc)
  expect_equal(nrow(res), 2L)
  expect_true(all(res$approach == "B"))
  expect_true(all(is.finite(res$RL)))
})

test_that("approach B requires mc_sample", {
  s <- make_stationary_fit()
  expect_error(
    marginalize(s$fit, s$df, approaches = "B"),
    regexp = "mc_sample"
  )
})

test_that("missing year column throws informative error", {
  s   <- make_stationary_fit()
  df2 <- s$df
  df2$year <- NULL
  expect_error(
    marginalize(s$fit, df2, approaches = "C"),
    regexp = "year"
  )
})

test_that("interactions are passed through correctly", {
  set.seed(3L)
  n  <- 500L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(40, 5, 20)),
    a    = stats::rnorm(n + 40L),
    b    = stats::rnorm(n + 40L),
    year = rep(1981:2020, each = 13L)[seq_len(n + 40L)]
  )
  df$ab <- df$a * df$b
  fit <- fit_nhpp(df, threshold = 4, loc_vars = c("a", "b", "ab"),
                  lambda = 0, verbose = FALSE)
  ints <- list(ab = c("a", "b"))
  sc   <- list(mid = list(a = 0, b = 0))
  res  <- marginalize(fit, df, TRs = 10, approaches = "A",
                      scenarios = sc, interactions = ints)
  expect_true(is.finite(res$RL))
})

test_that("approach C reconstructs active interactions from their inputs", {
  set.seed(31L)
  n_years <- 20L
  n_obs <- 13L
  n <- n_years * n_obs
  df <- data.frame(
    y = c(stats::rexp(n - 30L, 0.3), stats::runif(30L, 5, 20)),
    a = stats::rnorm(n),
    b = stats::rnorm(n),
    year = rep(seq_len(n_years), each = n_obs)
  )
  df$ab <- df$a * df$b
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "ab",
                  penalty = "none", obs_per_year = n_obs,
                  verbose = FALSE)
  fit$par["mu.ab"] <- 0.2
  fit$par_oper <- fit$par
  res <- marginalize(
    fit, df, TRs = 10, approaches = "C", n_obs = n_obs,
    n_boot = 20L, interactions = list(ab = c("a", "b")), seed = 1L
  )
  expect_true(is.finite(res$RL))
})

test_that("approach B rejects missing interaction inputs", {
  set.seed(32L)
  n <- 260L
  df <- data.frame(
    y = c(stats::rexp(n - 30L, 0.3), stats::runif(30L, 5, 20)),
    a = stats::rnorm(n), b = stats::rnorm(n),
    year = rep(seq_len(20L), each = 13L)
  )
  df$ab <- df$a * df$b
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "ab",
                  penalty = "none", obs_per_year = 13,
                  verbose = FALSE)
  fit$par["mu.ab"] <- 0.2
  fit$par_oper <- fit$par
  expect_error(
    marginalize(
      fit, df, TRs = 10, approaches = "B", n_obs = 13L,
      mc_sample = list(data.frame(a = rnorm(13L))),
      interactions = list(ab = c("a", "b"))
    ),
    regexp = "missing required covariates"
  )
})

test_that("marginalize rejects invalid inputs", {
  s <- make_stationary_fit()
  expect_error(marginalize(list(), s$df), regexp = "must be an nhpp_fit object")
  expect_error(marginalize(s$fit, list()), regexp = "must be a data frame")
  expect_error(marginalize(s$fit, s$df, TRs = 1), regexp = "greater than 1")
  expect_error(marginalize(s$fit, s$df, approaches = "D"), regexp = "only 'A', 'B', or 'C'")
  expect_error(
    marginalize(s$fit, s$df, approaches = "B", mc_sample = data.frame(x = 1)),
    regexp = "non-empty list"
  )
  expect_error(marginalize(s$fit, s$df, n_obs = 0), regexp = "n_obs")
  expect_error(marginalize(s$fit, s$df, n_boot = 0), regexp = "n_boot")
  expect_error(
    marginalize(s$fit, s$df, z_hi = s$fit$threshold),
    regexp = "greater than the threshold"
  )

  bad_y <- s$df
  bad_y$y <- NA_real_
  expect_error(marginalize(s$fit, bad_y), regexp = "data\\$y")
})

test_that("return-level search handles boundary and failed brackets", {
  at_threshold <- margEVT:::.find_return_level(
    TR = 2, f_annual = function(z) 0.5,
    z_lo = 10, z_hi = 20
  )
  expect_equal(at_threshold, 10)

  no_upper_crossing <- margEVT:::.find_return_level(
    TR = 10, f_annual = function(z) 0.1,
    z_lo = 10, z_hi = 20
  )
  expect_true(is.na(no_upper_crossing))

  failed_evaluation <- margEVT:::.find_return_level(
    TR = 10, f_annual = function(z) stop("failure"),
    z_lo = 10, z_hi = 20
  )
  expect_true(is.na(failed_evaluation))

  nonfinite_after_expansion <- margEVT:::.find_return_level(
    TR = 10,
    f_annual = function(z) if (z <= 20) 0.1 else NaN,
    z_lo = 10, z_hi = 20
  )
  expect_true(is.na(nonfinite_after_expansion))
})

test_that("approach C preserves the caller's random-number stream", {
  s <- make_stationary_fit()
  set.seed(99L)
  expected <- runif(2L)

  set.seed(99L)
  first <- runif(1L)
  marginalize(s$fit, s$df, TRs = 10, approaches = "C",
              n_boot = 10L, seed = 123L)
  second <- runif(1L)

  expect_equal(c(first, second), expected)
})

test_that("marginalize approach C handles under-observed years", {
  # MUST use the covariate model. The stationary model skips the resampling block!
  s <- make_cov_fit()

  # The dataset has exactly 13 rows per year.
  # By passing n_obs = 50L, nrow(df_yr) < 50 evaluates to TRUE for all years.
  # This returns NULL for all bootstraps and throws the expected error flawlessly.
  expect_error(
    marginalize(s$fit, s$df, TRs = 10, approaches = "C", n_boot = 5L, n_obs = 50L),
    regexp = "approach C produced no valid bootstrap years"
  )
})
