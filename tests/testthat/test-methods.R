# Tests for S3 methods: summary object, logLik/nobs/AIC/BIC, predict,
# residuals, path and Q-Q plots, and small exported helpers.

make_methods_fit <- function(penalty = "none", lambda = 0, seed = 21L) {
  set.seed(seed)
  n <- 3000L
  df <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  df$y <- 10 + 1.5 * df$x + stats::rexp(n, 1 / 3)
  list(df = df,
       fit = fit_nhpp(df, threshold = 17, loc_vars = c("x", "z"),
                      penalty = penalty, lambda = lambda, verbose = FALSE))
}

test_that("summary returns a summary.nhpp_fit object with a coefficient table", {
  s <- summary(make_methods_fit()$fit)
  expect_s3_class(s, "summary.nhpp_fit")
  expect_true(all(c("estimate", "operational", "active") %in%
                    names(s$coefficients)))
  expect_true(is.finite(s$bic))
  expect_output(print(s), "Exceedances")
  expect_output(print(s), "KKT violation")
})

test_that("logLik, nobs, AIC and BIC are consistent with bic_nhpp", {
  fit <- make_methods_fit()$fit
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -fit$nllh_oper)
  expect_equal(nobs(fit), fit$n_exc)
  expect_equal(BIC(fit), bic_nhpp(fit))
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
})

test_that("predict.nhpp_fit wraps predict_params", {
  fit <- make_methods_fit()$fit
  nd <- data.frame(x = c(-1, 1), z = c(0, 0))
  expect_identical(predict(fit, nd), predict_params(fit, nd))
  expect_identical(predict(fit), predict_params(fit))
})

test_that("interarrival residuals are approximately standard exponential", {
  fit <- make_methods_fit()$fit
  z <- residuals(fit)
  expect_length(z, fit$n_exc - 1L)
  expect_true(all(z > 0))
  expect_equal(mean(z), 1, tolerance = 0.2)
  expect_equal(stats::var(z), 1, tolerance = 0.35)
})

test_that("interarrival residuals reproduce the discretized definition", {
  fit <- make_methods_fit()$fit
  prm <- predict_params(fit)
  tau <- margEVT:::.tail_measure_at_level(fit$threshold, prm$mu, prm$sigma, prm$xi)
  exc <- which(fit$y > fit$threshold)
  manual <- vapply(2:3, function(k)
    sum(tau[exc[k - 1L]:(exc[k] - 1L)]) / fit$obs_per_year, numeric(1L))
  expect_equal(residuals(fit)[1:2], manual)
})

test_that("magnitude residuals handle general and Gumbel branches", {
  fit <- make_methods_fit()$fit
  r <- residuals(fit, type = "magnitude")
  expect_length(r, fit$n_exc)
  expect_true(all(r >= 0))
  fit_g <- fit
  fit_g$fitted_oper$xi[] <- 0
  exc <- which(fit$y > fit$threshold)
  beta <- fit_g$fitted_oper$sigma[exc] +
    0 * (fit$threshold - fit_g$fitted_oper$mu[exc])
  expect_equal(residuals(fit_g, type = "magnitude"),
               (fit$y[exc] - fit$threshold) / beta)
})

test_that("residuals validate stored series, support and sample size", {
  fit <- make_methods_fit()$fit
  legacy <- fit
  legacy$y <- NULL
  expect_error(residuals(legacy), "does not store the observed series")
  bad <- fit
  bad$fitted_oper$xi[] <- 1
  bad$fitted_oper$mu[] <- fit$threshold + 1000
  expect_error(residuals(bad), "outside the fitted support")
  one <- fit
  one$y <- c(fit$threshold + 1, rep(fit$threshold - 1, length(fit$y) - 1L))
  expect_length(residuals(one), 0L)
})

test_that("plot.nhpp_fit draws Q-Q and path diagnostics", {
  fit <- make_methods_fit()$fit
  expect_silent(plot(fit, type = "qq"))
  expect_silent(plot(fit, type = "qq", residuals = "magnitude"))
  expect_error(plot(fit, type = "path"), "no regularization path")
  few <- fit
  few$y <- c(fit$threshold + 1, rep(fit$threshold - 1, length(fit$y) - 1L))
  expect_error(plot(few, type = "qq"), "too few residuals")
  set.seed(3L)
  df <- make_methods_fit()$df
  fit_bic <- fit_nhpp(df, threshold = 17, loc_vars = c("x", "z"),
                      penalty = "lasso", lambda = "bic",
                      lambda_grid = c(20, 5, 1, 0.2), verbose = FALSE)
  expect_silent(plot(fit_bic, type = "path"))
})

test_that("n_exceedances uses the stored series when y is missing", {
  fit <- make_methods_fit()$fit
  expect_equal(n_exceedances(fit), fit$n_exc)
  legacy <- fit
  legacy$y <- NULL
  expect_error(n_exceedances(legacy), "must be supplied")
})

test_that("lambda formatting handles scalar, named and unnamed vectors", {
  expect_identical(margEVT:::.format_lambda(2), "2")
  expect_identical(margEVT:::.format_lambda(c(mu = 1, sigma = 2, xi = 3)),
                   "mu=1, sigma=2, xi=3")
  expect_identical(margEVT:::.format_lambda(c(1, 2, 3)), "mu=1, sigma=2, xi=3")
})

test_that("seed, core and harmonic helpers behave as documented", {
  expect_identical(margEVT:::.maybe_seeded(NULL, 5), 5)
  a <- margEVT:::.maybe_seeded(1L, stats::runif(1L))
  b <- margEVT:::.maybe_seeded(1L, stats::runif(1L))
  expect_identical(a, b)
  expect_error(margEVT:::.validate_seed("a", "f"), "seed")
  expect_null(margEVT:::.validate_seed(NULL, "f"))
  expect_error(margEVT:::.validate_cores(0, "f"), "cores")
  expect_identical(margEVT:::.validate_cores(2, "f"), 2L)
  expect_equal(margEVT:::.is_harmonic(c("cos1", "sin2", "sen3", "x", "cos")),
               c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that(".par_lapply gives identical results serially and in parallel", {
  f <- function(i) i^2
  expect_identical(margEVT:::.par_lapply(1:3, f, cores = 1L), as.list((1:3)^2))
  expect_identical(margEVT:::.par_lapply(1:3, f, cores = 2L), as.list((1:3)^2))
})
