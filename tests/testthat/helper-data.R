# Simulated data shared by several test files (testthat sources helper-*.R
# files before running the tests).

# Thirty years of 50 observations: a location covariate `x`, an irrelevant
# covariate `z` and a trend normalized to [0, 1].
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

# Autocorrelated location covariate, 24 years of 50 observations, and an
# unpenalized fit.
make_bt_df <- function(seed = 3L, n_years = 24L, n_obs = 50L) {
  set.seed(seed)
  n <- n_years * n_obs
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), n))
  df <- data.frame(x = x, year = rep(seq_len(n_years) + 1990L, each = n_obs))
  df$y <- 5 + 0.8 * df$x + stats::rexp(n, 1 / 2)
  fit <- fit_nhpp(df, threshold = 9, loc_vars = "x", penalty = "none",
                  obs_per_year = n_obs, verbose = FALSE)
  list(df = df, fit = fit, n_obs = n_obs)
}

# One relevant (`x1`) and two irrelevant covariates, for penalized estimation.
make_pen_df <- function(n = 1500L, seed = 4L) {
  set.seed(seed)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   x3 = stats::rnorm(n))
  df$y <- 10 + 2 * df$x1 + stats::rexp(n, 1 / 3)
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
