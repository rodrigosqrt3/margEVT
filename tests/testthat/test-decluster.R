# Tests for the extremal index and declustering.

make_ar_series <- function(seed = 3L, n = 5000L, ar = 0.7) {
  set.seed(seed)
  y <- as.numeric(stats::arima.sim(list(ar = ar), n))
  list(y = y, u = unname(stats::quantile(y, 0.95)))
}

test_that("runs declustering keeps one maximum per cluster", {
  y <- c(0, 5, 7, 0, 0, 6, 0, 8, 8, 0, 0, 0, 9)
  d1 <- decluster(y, 4, method = "runs", run_length = 1L)
  expect_identical(attr(d1, "clusters"), c(1L, 1L, 2L, 3L, 3L, 4L))
  expect_identical(as.numeric(d1), c(0, 4, 7, 0, 0, 6, 0, 8, 4, 0, 0, 0, 9))
  d2 <- decluster(y, 4, method = "runs", run_length = 2L)
  expect_equal(attr(d2, "n_clusters"), 3)
  expect_identical(which(d2 > 4), c(3L, 8L, 13L))
  d0 <- decluster(y, 4, method = "runs", run_length = 0L, replace_with = -1)
  expect_identical(as.numeric(d0), y)
})

test_that("extremal_index and decluster agree with each other", {
  s <- make_ar_series()
  ei <- extremal_index(s$y, s$u)
  expect_named(ei, c("theta", "n_clusters", "run_length"))
  expect_gt(ei[["theta"]], 0)
  expect_lt(ei[["theta"]], 1)
  yd <- decluster(s$y, s$u)
  expect_identical(sum(yd > s$u), as.integer(ei[["n_clusters"]]))
  expect_equal(attr(yd, "run_length"), ei[["run_length"]])
  er <- extremal_index(s$y, s$u, method = "runs", run_length = 3L)
  expect_equal(er[["theta"]], er[["n_clusters"]] / sum(s$y > s$u))
  # independent series: theta close to 1, small run length
  set.seed(8L)
  z <- stats::rnorm(20000L)
  ez <- extremal_index(z, stats::quantile(z, 0.98))
  expect_gt(ez[["theta"]], 0.85)
})

test_that("results match extRemes", {
  skip_if_not_installed("extRemes")
  s <- make_ar_series()
  ex <- extRemes::extremalindex(s$y, s$u, method = "intervals")
  ei <- extremal_index(s$y, s$u)
  expect_equal(ei[["theta"]], unname(ex[1L]))
  expect_equal(ei[["run_length"]], unname(ex[3L]))
  dr <- extRemes::decluster(s$y, s$u, method = "runs", r = ex[3L])
  expect_equal(as.numeric(decluster(s$y, s$u)), as.numeric(c(dr)))
  d1 <- extRemes::decluster(s$y, s$u, method = "runs", r = 1)
  expect_equal(as.numeric(decluster(s$y, s$u, method = "runs")),
               as.numeric(c(d1)))
})

test_that("decluster works as the backtest decluster function", {
  s <- make_ar_series()
  f <- function(v) decluster(v, s$u)
  out <- f(s$y)
  expect_true(is.numeric(out) && length(out) == length(s$y))
})

test_that("intervals run length covers the edge cases", {
  expect_identical(margEVT:::.run_length_intervals(c(1, 1, 1), 1), 0L)
  expect_identical(margEVT:::.run_length_intervals(c(5, 1, 9, 1), 0.4), 5L)
  expect_identical(margEVT:::.run_length_intervals(c(5, 1, 9, 1), 0.5), 1L)
  expect_equal(margEVT:::.theta_intervals(c(1, 1, 1)), 1)
})

test_that("decluster and extremal_index validate their input", {
  expect_error(decluster(c(1, 2, NA), 0), "finite")
  expect_error(decluster(c(1, 2, 3), c(0, 1)), "threshold")
  expect_error(decluster(c(1, 2, 3), 2.5), "two exceedances")
  expect_error(decluster(c(1, 5, 6), 2, replace_with = 3), "replace_with")
  expect_error(extremal_index(c(1, 5, 6), 2, method = "runs",
                              run_length = -1), "run_length")
})
