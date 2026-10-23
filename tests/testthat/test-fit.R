make_test_df <- function(n = 300L, seed = 42L) {
  set.seed(seed)
  data.frame(
    y = c(stats::rexp(n, rate = 0.5), stats::runif(30L, 5, 15)),
    x = stats::rnorm(n + 30L)
  )
}

test_that("fit_nhpp returns an nhpp_fit object", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_s3_class(fit, "nhpp_fit")
  expect_equal(fit$active_tol, 1e-2)
})

test_that("fit_nhpp stores a user-supplied activity tolerance", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none",
                  active_tol = 5e-3, verbose = FALSE)
  expect_equal(fit$active_tol, 5e-3)
})

test_that("fit_nhpp stores smooth and operational coefficient vectors", {
  df <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 500,
                  active_tol = 1e-2, verbose = FALSE)
  expect_named(fit$par_oper, names(fit$par))
  expected_x <- if (abs(fit$par["mu.x"]) < fit$active_tol) 0 else
    unname(fit$par["mu.x"])
  expect_equal(unname(fit$par_oper["mu.x"]), expected_x)
  expect_true(is.finite(fit$nllh_oper))
  expect_named(coef(fit, operational = TRUE), names(fit$par))
})

test_that("penalty=none gives lambda=0 and correct parameter names", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$lambda, 0)
  expect_named(fit$par,
               c("mu.(Intercept)", "sigma.(Intercept)", "xi.(Intercept)"))
})

test_that("penalty=lasso forces alpha=1", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 1)
  expect_equal(fit$penalty, "lasso")
})

test_that("penalty=ridge forces alpha=0", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "ridge", lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 0)
  expect_equal(fit$penalty, "ridge")
})

test_that("penalty=elnet uses supplied alpha", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "elnet", alpha = 0.3,
                  lambda = 1, verbose = FALSE)
  expect_equal(fit$alpha, 0.3)
})

test_that("large lambda shrinks covariate to zero under lasso", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 500,
                  verbose = FALSE)
  expect_lt(abs(fit$par["mu.x"]), 0.01)
})

test_that("free_vars not shrunk by large lambda", {
  df       <- make_test_df()
  df$trend <- seq(0, 1, length.out = nrow(df))
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars  = c("trend", "x"),
                  free_vars = "trend",
                  penalty   = "lasso", lambda = 500,
                  verbose   = FALSE)
  expect_lt(abs(fit$par["mu.x"]),    0.01)
  expect_gt(abs(fit$par["mu.trend"]), 0)
})

test_that("numeric lambda skips grid search", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = "x",
                  penalty = "lasso", lambda = 0.5,
                  verbose = FALSE)
  expect_equal(fit$lambda, 0.5)
})

test_that("coef() returns parameter vector", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_identical(coef(fit), fit$par)
})

test_that("fitted values have correct length and sigma positive", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
  expect_length(fit$fitted$mu,    nrow(df))
  expect_length(fit$fitted$sigma, nrow(df))
  expect_length(fit$fitted$xi,    nrow(df))
  expect_true(all(fit$fitted$sigma > 0))
})

test_that("calc_hessian stores a square matrix", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none",
                  calc_hessian = TRUE, verbose = FALSE)
  expect_true(is.matrix(fit$hessian))
  expect_equal(nrow(fit$hessian), length(fit$par))
})

test_that("fit_nhpp input validation throws expected errors", {
  df <- make_test_df()
  expect_error(fit_nhpp(list(y = 1:10), threshold = 4), regexp = "must be a data frame")
  expect_error(fit_nhpp(data.frame(x = 1:10), threshold = 4), regexp = "must contain a column named `y`")
  expect_error(fit_nhpp(df, threshold = c(1, 4)), regexp = "single numeric value")
  expect_error(fit_nhpp(df, threshold = 4, active_tol = 0), regexp = "active_tol")
  expect_error(fit_nhpp(df, threshold = 4, penalty = "lasso", lambda = "invalid"), regexp = "numeric value or")
  expect_error(fit_nhpp(df, threshold = 4, penalty = "lasso", lambda = -1), regexp = "lambda")
  expect_error(
    fit_nhpp(df, threshold = 4, penalty = "elnet", alpha = 2, lambda = 1),
    regexp = "alpha"
  )
  expect_error(
    fit_nhpp(df, threshold = 4, penalty = "elnet",
             alpha = c(0.5, 0.5, 0.5), lambda = 1),
    regexp = "must be named"
  )
  df_na <- df
  df_na$x[1L] <- NA_real_
  expect_error(
    fit_nhpp(df_na, threshold = 4, loc_vars = "x", penalty = "none"),
    regexp = "covariates"
  )
  expect_error(fit_nhpp(df, threshold = 4, maxit = 0), regexp = "maxit")
  expect_error(
    fit_nhpp(df, threshold = 4, obs_per_year = 0),
    regexp = "obs_per_year"
  )
  expect_error(
    fit_nhpp(df, threshold = 4, penalize_shape = NA),
    regexp = "penalize_shape"
  )
  expect_error(
    fit_nhpp(df, threshold = 4, calc_hessian = NA),
    regexp = "calc_hessian"
  )
  expect_error(
    fit_nhpp(df, threshold = 4, verbose = NA),
    regexp = "verbose"
  )

  df_bad_y <- df
  df_bad_y$y[1L] <- NA_real_
  expect_error(fit_nhpp(df_bad_y, threshold = 4), regexp = "df\\$y")

  df_no_exc <- df
  df_no_exc$y <- 0
  expect_error(fit_nhpp(df_no_exc, threshold = 4), regexp = "no observations exceed")

  expect_error(
    fit_nhpp(
      df, threshold = 4, loc_vars = "x", penalty = "elnet",
      alpha = c(a = 0.5, b = 0.5, c = 0.5), lambda = 1
    ),
    regexp = "named mu, sigma, and xi"
  )
})

test_that("fit_nhpp accepts named block-specific controls", {
  df <- make_test_df()
  fit <- fit_nhpp(
    df, threshold = 4, loc_vars = "x", penalty = "elnet",
    alpha = c(xi = 0.3, mu = 0.5, sigma = 0.4),
    lambda = c(sigma = 0.2, xi = 0.1, mu = 0.3),
    verbose = FALSE
  )
  expect_equal(names(fit$alpha), c("mu", "sigma", "xi"))
  expect_equal(names(fit$lambda), c("mu", "sigma", "xi"))
})

test_that("common BIC scaling applies the same lambda to all blocks", {
  df <- make_test_df()
  fit <- fit_nhpp(
    df, threshold = 4, loc_vars = "x", scale_vars = "x",
    penalty = "lasso", lambda = "bic", lambda_scaling = "common",
    verbose = FALSE
  )
  expect_equal(unname(fit$lambda), rep(unname(fit$lambda[1L]), 3L))
  expect_identical(fit$lambda_scaling, "common")
})

test_that("fit_nhpp validates lambda scaling", {
  df <- make_test_df()
  expect_error(
    fit_nhpp(df, threshold = 4, lambda_scaling = "unknown"),
    regexp = "arg"
  )
})

test_that(".fit_at_lambda returns a failed result when both optimizers error", {
  df <- make_test_df()
  dm <- margEVT:::build_design_matrices(
    df, loc_vars = "x", scale_vars = NULL,
    shape_vars = NULL, free_vars = NULL
  )
  attr(dm, "threshold") <- 4
  init <- rep(0, ncol(dm$X_mu) + ncol(dm$X_sigma) + ncol(dm$X_xi))

  testthat::local_mocked_bindings(
    .optim_nhpp = function(...) stop("optimizer failure")
  )

  result <- margEVT:::.fit_at_lambda(
    dm = dm, y = df$y, threshold = 4,
    lambda = 0, alpha = 1, penalize_shape = TRUE,
    init = init
  )

  expect_false(result$converged)
  expect_identical(result$par, init)
  expect_true(is.na(result$nllh_pen))
  expect_true(is.na(result$nllh_raw))
  expect_null(result$hessian)
})

test_that("fit_nhpp throws warning on < 5 exceedances", {
  df_few <- data.frame(y = c(1, 2, 3, 4.5, 4.5, 4.5), x = rnorm(6))
  # With three exceedances the fit is degenerate; only the small-sample
  # warning is under test here, so the convergence warning is muffled.
  expect_warning(
    withCallingHandlers(
      fit_nhpp(df_few, threshold = 4, penalty = "none", verbose = TRUE),
      warning = function(w)
        if (grepl("did not converge", conditionMessage(w)))
          invokeRestart("muffleWarning")),
    regexp = "fewer than 5 exceedances"
  )
})

test_that("fit_nhpp throws warning on non-convergence", {
  df <- make_test_df()
  expect_warning(
    fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none", maxit = 1L, verbose = TRUE),
    regexp = "optimizer did not converge"
  )
})

# -----------------------------------------------------------------------------
# S3 Print Method Test
# -----------------------------------------------------------------------------
test_that("print.nhpp_fit prints without error", {
  df  <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)

  expect_output(print(fit), regexp = "-- nhpp_fit --")
  expect_output(print(fit), regexp = "Threshold")
  expect_output(print(fit), regexp = "Active params")
  expect_output(print(fit), regexp = "Coefficients")
})

test_that(".warmstart handles empty pen_idx_all", {
  df <- make_test_df()
  fit <- fit_nhpp(df, threshold = 4, loc_vars = NULL, scale_vars = NULL,
                  penalize_shape = FALSE, penalty = "lasso", lambda = 1.0, verbose = FALSE)
  expect_s3_class(fit, "nhpp_fit")
})

test_that(".fit_at_lambda handles missing init correctly", {
  df <- make_test_df()
  dm <- margEVT:::build_design_matrices(df, loc_vars = "x", scale_vars = NULL, shape_vars = NULL, free_vars = NULL)
  attr(dm, "threshold") <- 4

  res <- margEVT:::.fit_at_lambda(dm, df$y, threshold = 4, lambda = 0, alpha = 1,
                                  penalize_shape = TRUE, init = NULL)

  expect_type(res, "list")
  expect_true(res$converged)
})

test_that(".warmstart handles all(lambda == 0) early exit", {
  init_mock <- c(1.0, 2.0, 3.0)
  res <- margEVT:::.warmstart(init = init_mock, dm = NULL, y = NULL,
                              threshold = 0, lambda = 0, alpha = 1,
                              penalize_shape = TRUE)
  expect_identical(res, init_mock)
})

test_that(".warmstart applies coordinate-specific soft thresholds", {
  df <- make_pen_df(500L)
  dm <- build_design_matrices(df, loc_vars = "x2", scale_vars = "x3")
  init <- c(15, 0.4, log(3), 0.3, 0.05)
  ws <- margEVT:::.warmstart(init, dm, df$y, 16,
                             lambda = c(mu = 1e6, sigma = 0, xi = 0),
                             alpha = 1, penalize_shape = TRUE)
  expect_equal(ws[2L], 0)
  expect_equal(ws[4L], init[4L])
  expect_identical(
    margEVT:::.warmstart(init, dm, df$y, 16, lambda = c(mu = 1, sigma = 1, xi = 1),
                         alpha = 1, penalize_shape = TRUE,
                         penalty_weights = rep(0, 5L)),
    init)
})

# -----------------------------------------------------------------------------
# Optimizer fallbacks and KKT certification
# -----------------------------------------------------------------------------

test_that(".fit_at_lambda tries a BFGS restart and L-BFGS-B when KKT fails", {
  df <- make_pen_df(300L)
  dm <- build_design_matrices(df, loc_vars = "x1")
  calls <- character(0L)
  testthat::local_mocked_bindings(
    .optim_nhpp = function(par, fn, gr, method, control) {
      calls <<- c(calls, method)
      list(par = par, value = fn(par), convergence = 0L)
    }
  )
  res <- margEVT:::.fit_at_lambda(dm, df$y, 16, lambda = 0, alpha = 1,
                                  penalize_shape = TRUE,
                                  init = c(15, 0, log(3), 0.05))
  expect_identical(calls, c("BFGS", "BFGS", "L-BFGS-B"))
  expect_false(res$converged)
  expect_true(is.finite(res$kkt))
})

test_that(".fit_at_lambda discards optimizer results outside the support", {
  df <- make_pen_df(300L)
  dm <- build_design_matrices(df)
  testthat::local_mocked_bindings(
    .optim_nhpp = function(par, fn, gr, method, control)
      list(par = par, value = 1e9, convergence = 0L)
  )
  res <- margEVT:::.fit_at_lambda(dm, df$y, 16, lambda = 0, alpha = 1,
                                  penalize_shape = TRUE,
                                  init = c(15, log(3), 0.05))
  expect_false(res$converged)
  expect_true(is.na(res$kkt))
})

# -----------------------------------------------------------------------------
# Penalty weights and starting values
# -----------------------------------------------------------------------------

test_that("adaptive weights are stored and applied", {
  df <- make_pen_df()
  expect_message(
    fit <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2", "x3"),
                    penalty = "lasso", lambda = 2,
                    penalty_weights = "adaptive", verbose = TRUE),
    "adaptive weights from mle")
  w <- fit$penalty_weights
  expect_named(w, names(fit$par))
  expect_equal(unname(w[c("mu.(Intercept)", "sigma.(Intercept)",
                          "xi.(Intercept)")]), c(0, 0, 0))
  expect_true(all(w[c("mu.x1", "mu.x2", "mu.x3")] > 0))
  expect_lt(w["mu.x1"], w["mu.x2"])
  expect_identical(fit$weights_source, "adaptive: mle")
  expect_output(print(fit), "Weights")
  expect_output(print(summary(fit)), "Weights")
})

test_that("adaptive weights fall back to ridge and fail informatively", {
  df <- make_pen_df(400L)
  dm <- build_design_matrices(df, loc_vars = "x1")
  init <- c(15, 0, log(3), 0.05)
  testthat::local_mocked_bindings(
    .fit_at_lambda = function(dm, y, threshold, lambda, alpha, ...)
      list(converged = lambda >= 1, par = c(15, 0.5, log(3), 0.05))
  )
  aw <- margEVT:::.adaptive_weights(dm, df$y, 16, init, TRUE, 365.25)
  expect_identical(aw$source, "ridge(lambda=1)")
  expect_equal(aw$weights[2L], 2)
  expect_error(
    margEVT:::.adaptive_weights(dm, df$y, 16, init, TRUE, 365.25,
                                ridge_lambdas = 0.5),
    "convergent initial fit")
})

test_that("numeric penalty weights and init are validated and reordered", {
  df <- make_pen_df(600L)
  base <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2"),
                   penalty = "none", verbose = FALSE)
  nm <- names(base$par)
  w <- stats::setNames(c(0, 1, 5, 0, 0), nm)
  fit <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2"),
                  penalty = "lasso", lambda = 1, penalty_weights = rev(w),
                  init = rev(base$par), verbose = FALSE)
  expect_identical(fit$weights_source, "user")
  expect_equal(fit$penalty_weights, w)
  expect_error(fit_nhpp(df, 16, loc_vars = "x1", penalty = "lasso", lambda = 1,
                        penalty_weights = c(1, -1, 0, 0), verbose = FALSE),
               "non-negative")
  expect_error(fit_nhpp(df, 16, loc_vars = "x1", penalty = "lasso", lambda = 1,
                        penalty_weights = c(a = 1, b = 1, c = 1, d = 1),
                        verbose = FALSE), "names of `penalty_weights`")
  expect_error(fit_nhpp(df, 16, loc_vars = "x1", penalty = "lasso", lambda = 1,
                        penalty_weights = "other", verbose = FALSE),
               "penalty_weights")
  expect_error(fit_nhpp(df, 16, loc_vars = "x1", init = 1:2, verbose = FALSE),
               "init")
  expect_error(fit_nhpp(df, 16, loc_vars = "x1",
                        init = c(a = 1, b = 1, c = 1, d = 1), verbose = FALSE),
               "names of `init`")
  expect_error(fit_nhpp(df, 16, grad_tol = 0, verbose = FALSE), "grad_tol")
  expect_error(fit_nhpp(df, 16, adaptive_gamma = -1, verbose = FALSE),
               "adaptive_gamma")
})

test_that("adaptive weights are ignored for unpenalized fits", {
  df <- make_pen_df(400L)
  fit <- fit_nhpp(df, threshold = 16, loc_vars = "x1", penalty = "none",
                  penalty_weights = "adaptive", verbose = FALSE)
  expect_null(fit$penalty_weights)
  expect_identical(fit$weights_source, "none")
})

# -----------------------------------------------------------------------------
# Lower bound on the shape parameter
# -----------------------------------------------------------------------------

test_that("fit_nhpp converges at the shape bound for a short-tailed sample", {
  set.seed(5)
  n <- 4000L
  df <- data.frame(y = stats::runif(n, 0, 10), t = seq(0, 1, length.out = n))
  # uniform tail: GPD shape -1, so the bound -0.5 is active
  fit <- fit_nhpp(df, threshold = 8, penalty = "none", verbose = FALSE)
  expect_true(fit$converged)
  expect_true(fit$shape_at_bound)
  expect_equal(unname(fit$par["xi.(Intercept)"]), -0.5, tolerance = 1e-6)
  expect_output(print(fit), "Shape bound  : active \\(xi >= -0.5\\)")
  fit_t <- fit_nhpp(df, threshold = 8, shape_vars = "t", penalty = "lasso",
                    lambda = 1, verbose = FALSE)
  expect_true(fit_t$converged)
  expect_gte(min(fit_t$fitted$xi), -0.5 - 1e-8)
  expect_error(fit_nhpp(df, threshold = 8, shape_lower = 0.1, verbose = FALSE),
               "shape_lower")
})

test_that("fits away from the bound are unaffected by it", {
  set.seed(6)
  df <- data.frame(y = stats::rexp(3000L, 1 / 3))
  f1 <- fit_nhpp(df, threshold = 8, penalty = "none", verbose = FALSE)
  f2 <- fit_nhpp(df, threshold = 8, penalty = "none", shape_lower = -Inf,
                 verbose = FALSE)
  expect_false(f1$shape_at_bound)
  expect_equal(f1$par, f2$par, tolerance = 1e-8)
})

test_that("refits that start at the shape bound converge (active set)", {
  set.seed(5)
  n <- 4000L
  df <- data.frame(y = stats::runif(n, 0, 10), t = seq(0, 1, length.out = n))
  fit <- fit_nhpp(df, threshold = 8, shape_vars = "t", penalty = "lasso",
                  lambda = 1, verbose = FALSE)
  expect_true(fit$shape_at_bound)
  df$y <- stats::runif(n, 0, 10)                        # new sample, same law
  calls <- character(0)
  testthat::local_mocked_bindings(
    .optim_nhpp = function(...) {
      a <- list(...); calls <<- c(calls, a$method); stats::optim(...)
    }
  )
  res <- margEVT:::.fit_at_lambda(fit$dm, df$y, 8, lambda = fit$lambda,
                                  alpha = 1, penalize_shape = TRUE,
                                  init = fit$par, obs_per_year = 365.25,
                                  shape_lower = -0.5)
  expect_true(res$converged)
  expect_identical(calls[1L], "BFGS")                   # active set, not L-BFGS-B
  expect_gte(min(fit$dm$X_xi %*% utils::tail(res$par, 2L)), -0.5 - 1e-8)
})

test_that("active-set failure at the bound falls back to bounded L-BFGS-B", {
  set.seed(5)
  n <- 4000L
  df <- data.frame(y = stats::runif(n, 0, 10), t = seq(0, 1, length.out = n))
  fit <- fit_nhpp(df, threshold = 8, shape_vars = "t", penalty = "lasso",
                  lambda = 1, verbose = FALSE)
  expect_true(fit$shape_at_bound)
  calls <- character(0)
  testthat::local_mocked_bindings(
    .optim_nhpp = function(...) {
      a <- list(...); calls <<- c(calls, a$method)
      if (identical(a$method, "BFGS")) stop("forced BFGS failure")
      stats::optim(...)
    }
  )
  res <- margEVT:::.fit_at_lambda(fit$dm, df$y, 8, lambda = fit$lambda,
                                  alpha = 1, penalize_shape = TRUE,
                                  init = fit$par, obs_per_year = 365.25,
                                  shape_lower = -0.5)
  expect_identical(calls[1L], "BFGS")                   # active set tried first
  expect_true("L-BFGS-B" %in% calls)
  expect_true(all(is.finite(res$par)))
  expect_gte(min(fit$dm$X_xi %*% utils::tail(res$par, 2L)), -0.5 - 1e-6)
})

# -----------------------------------------------------------------------------
# Parameter scaling for the optimizer
# -----------------------------------------------------------------------------

test_that("parameter scaling reaches a fit at least as good as unscaled", {
  set.seed(5L)
  n <- 3000L
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n), x3 = stats::rnorm(n))
  df$y <- 30 + 3 * df$x1 + exp(1.5 + 0.2 * df$x2) *
    (((100 * -log(stats::runif(n)))^0.1 - 1) / -0.1)
  u <- unname(stats::quantile(df$y, 0.95))
  args <- list(df = df, threshold = u, loc_vars = c("x1", "x2", "x3"),
               scale_vars = c("x1", "x2", "x3"), penalty = "lasso",
               lambda = 1, penalty_weights = "adaptive", obs_per_year = 100,
               verbose = FALSE)
  fa <- do.call(fit_nhpp, c(args, scaling = "auto"))
  fn <- do.call(fit_nhpp, c(args, scaling = "none"))
  expect_identical(fa$scaling, "auto")
  expect_identical(fn$scaling, "none")
  expect_true(fa$converged && fn$converged)
  expect_lte(fa$nllh_pen, fn$nllh_pen + 1e-6 * abs(fn$nllh_pen))
  expect_identical(fa$par_oper != 0, fn$par_oper != 0)
  expect_equal(unname(fa$par), unname(fn$par), tolerance = 1e-2)
  ps <- margEVT:::.par_scale(fa$par, fa$dm, df$y, u, 100)
  expect_true(all(is.finite(ps)) && all(ps > 0))
})

test_that(".par_scale returns NULL when the curvature cannot be computed", {
  testthat::local_mocked_bindings(
    .nhpp_eval = function(...) list(value = 1e9, grad = c(0, 0))
  )
  expect_null(margEVT:::.par_scale(c(0, 0), NULL, NULL, 0, 1))
})

test_that(".par_scale returns NULL when a perturbed point leaves the support", {
  n_calls <- 0L
  testthat::local_mocked_bindings(
    .nhpp_eval = function(...) {
      n_calls <<- n_calls + 1L
      if (n_calls == 1L) list(value = 1, grad = c(1, 1)) else
        list(value = 1e9, grad = c(NA_real_, NA_real_))
    }
  )
  expect_null(margEVT:::.par_scale(c(0, 0), NULL, NULL, 0, 1))
  expect_equal(n_calls, 3L)                             # e0 + one per coordinate
})
