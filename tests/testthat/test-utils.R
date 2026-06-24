# -----------------------------------------------------------------------------
# Test Helper
# -----------------------------------------------------------------------------
make_utils_fit <- function(seed = 1L) {
  set.seed(seed)
  n  <- 300L
  df <- data.frame(
    y = c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20)),
    x = stats::rnorm(n + 30L)
  )
  fit_nhpp(df, threshold = 4, loc_vars = "x",
           penalty = "none", verbose = FALSE)
}

# -----------------------------------------------------------------------------
# summary.nhpp_fit Tests
# -----------------------------------------------------------------------------
test_that("summary.nhpp_fit runs without error", {
  fit <- make_utils_fit()
  expect_output(summary(fit), regexp = "nhpp_fit summary")
})

test_that("summary.nhpp_fit shows threshold", {
  fit <- make_utils_fit()
  expect_output(summary(fit), regexp = "Threshold")
})

test_that("summary.nhpp_fit handles Hessian diagnostic branches (Mocked)", {
  # 1. Helper mock base list
  base_mock <- list(
    threshold = 4, penalty = "lasso", lambda = 0.1, alpha = 1,
    penalize_shape = TRUE, obs_per_year = 365.25, converged = TRUE,
    nllh_raw = 100, nllh_pen = 105, par = c(a = 1.0, b = 0.0) # b = 0 triggers the 'Shrunk to zero' print
  )
  class(base_mock) <- "nhpp_fit"

  # A. Test: Shrunk parameters printing & Hessian positive-definite
  fit_ok <- base_mock
  fit_ok$hessian <- matrix(c(2, 0, 0, 2), nrow = 2) # positive definite
  expect_output(summary(fit_ok), regexp = "Hessian positive definite")
  expect_output(summary(fit_ok), regexp = "Shrunk to zero : b")

  # B. Test: Hessian not positive-definite
  fit_npd <- base_mock
  fit_npd$hessian <- matrix(c(-1, 0, 0, 2), nrow = 2) # has negative eigenvalue
  expect_output(summary(fit_npd), regexp = "Hessian not positive definite")

  # C. Test: Hessian ill-conditioned
  fit_ill <- base_mock
  fit_ill$hessian <- matrix(c(1e9, 0, 0, 1e-3), nrow = 2) # extremely high condition number
  expect_output(summary(fit_ill), regexp = "Hessian ill-conditioned")
})

# -----------------------------------------------------------------------------
# is_nhpp_fit Tests
# -----------------------------------------------------------------------------
test_that("is_nhpp_fit returns TRUE for nhpp_fit", {
  fit <- make_utils_fit()
  expect_true(is_nhpp_fit(fit))
})

test_that("is_nhpp_fit returns FALSE for plain list", {
  expect_false(is_nhpp_fit(list(a = 1)))
})

# -----------------------------------------------------------------------------
# bic_nhpp Tests
# -----------------------------------------------------------------------------
test_that("bic_nhpp returns a finite scalar", {
  fit <- make_utils_fit()
  b   <- bic_nhpp(fit)
  expect_length(b, 1L)
  expect_true(is.finite(b))
})

test_that("bic_nhpp rejects non-nhpp_fit objects", {
  expect_error(bic_nhpp(list(a = 1)), regexp = "must be an nhpp_fit object")
})

test_that("bic_nhpp handles non-finite raw likelihood", {
  mock_inf <- list(nllh_raw = NA_real_)
  class(mock_inf) <- "nhpp_fit"
  expect_identical(bic_nhpp(mock_inf), NA_real_)
})

# -----------------------------------------------------------------------------
# n_exceedances Tests
# -----------------------------------------------------------------------------
test_that("n_exceedances returns correct count", {
  set.seed(1L)
  n  <- 300L
  y  <- c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20))
  df <- data.frame(y = y, x = stats::rnorm(n + 30L))
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_equal(n_exceedances(fit, y), sum(y > 4))
})

test_that("n_exceedances rejects non-nhpp_fit objects", {
  expect_error(n_exceedances(list(a = 1), 1:10), regexp = "must be an nhpp_fit object")
})

# -----------------------------------------------------------------------------
# rl_table Tests
# -----------------------------------------------------------------------------
test_that("rl_table produces wide format", {
  set.seed(1L)
  n  <- 300L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20)),
    x    = stats::rnorm(n + 30L),
    year = rep(1991:2020, each = 11L)[seq_len(n + 30L)]
  )
  fit  <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  marg <- marginalize(fit, df, TRs = c(10, 50), approaches = "A")
  wide <- rl_table(marg)
  expect_true(all(c("approach", "scenario", "T10", "T50") %in% names(wide)))
  expect_equal(nrow(wide), 1L)
})

test_that("rl_table with multiple scenarios has correct rows", {
  set.seed(1L)
  n  <- 300L
  df <- data.frame(
    y    = c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20)),
    x    = stats::rnorm(n + 30L),
    year = rep(1991:2020, each = 11L)[seq_len(n + 30L)]
  )
  fit  <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                   penalty = "none", verbose = FALSE)
  sc   <- list(low = list(x = -1), mid = list(x = 0), high = list(x = 1))
  marg <- marginalize(fit, df, TRs = c(10, 50),
                      approaches = "A", scenarios = sc)
  wide <- rl_table(marg)
  expect_equal(nrow(wide), 3L)
  expect_equal(sort(wide$scenario), c("high", "low", "mid"))
})

test_that("rl_table rejects wrong input", {
  expect_error(rl_table(list()), regexp = "data frame")
  expect_error(rl_table(data.frame(a = 1)), regexp = "columns")
})
