library(testthat)
library(margEVT)

# Helper: fit stationary model with synthetic data
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

# Helper: fit with one covariate
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

test_that("marginalize rejects invalid inputs", {
  s <- make_stationary_fit()
  expect_error(marginalize(list(), s$df), regexp = "must be an nhpp_fit object")
  expect_error(marginalize(s$fit, list()), regexp = "must be a data frame")
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
