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

test_that("summary.nhpp_fit runs without error", {
  fit <- make_utils_fit()
  expect_output(summary(fit), regexp = "nhpp_fit summary")
})

test_that("summary.nhpp_fit shows threshold", {
  fit <- make_utils_fit()
  expect_output(summary(fit), regexp = "Threshold")
})

test_that("is_nhpp_fit returns TRUE for nhpp_fit", {
  fit <- make_utils_fit()
  expect_true(is_nhpp_fit(fit))
})

test_that("is_nhpp_fit returns FALSE for plain list", {
  expect_false(is_nhpp_fit(list(a = 1)))
})

test_that("n_exceedances returns correct count", {
  set.seed(1L)
  n  <- 300L
  y  <- c(stats::rexp(n, 0.3), stats::runif(30L, 5, 20))
  df <- data.frame(y = y, x = stats::rnorm(n + 30L))
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_equal(n_exceedances(fit, y), sum(y > 4))
})

test_that("bic_nhpp returns a finite scalar", {
  fit <- make_utils_fit()
  b   <- bic_nhpp(fit)
  expect_length(b, 1L)
  expect_true(is.finite(b))
})

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
