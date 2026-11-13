# Tests for the proximal Newton solver (newton.R): analytical Hessian,
# coordinate descent, agreement with the BFGS solver, the shape bound and
# the fallbacks.

test_that("analytical Hessian matches finite differences of the gradient", {
  df <- make_pen_df()
  u <- 16
  fit <- fit_nhpp(df, threshold = u, loc_vars = c("x1", "x2"),
                  scale_vars = "x3", penalty = "none", obs_per_year = 50,
                  verbose = FALSE)
  dm <- fit$dm
  for (dxi in c(0, 0.1, -0.05)) {
    p <- unname(fit$par)
    p[6] <- p[6] + dxi
    h <- margEVT:::.nhpp_hessian(p, dm, df$y, u, 50, -0.5)
    num <- sapply(seq_along(p), function(k) {
      e <- 1e-5 * max(1, abs(p[k]))
      a <- p; b <- p; a[k] <- a[k] + e; b[k] <- b[k] - e
      (margEVT:::.nhpp_eval(a, dm, df$y, u, 50, TRUE)$grad -
         margEVT:::.nhpp_eval(b, dm, df$y, u, 50, TRUE)$grad) / (2 * e)
    })
    expect_equal(h$hessian, num, tolerance = 1e-5)
    # absolute comparison: at the MLE the gradient is close to zero, so a
    # relative tolerance would only measure rounding (it differs by platform)
    g <- margEVT:::.nhpp_eval(p, dm, df$y, u, 50, TRUE)$grad
    expect_lt(max(abs(h$grad - g)), 1e-8 * max(1, abs(g)))
  }
  expect_null(margEVT:::.nhpp_hessian(c(11, 0, 0, 0, 0, 1e-9), dm, df$y, u, 50,
                                      -0.5))                       # Gumbel
})

test_that("coordinate descent solves the penalized quadratic exactly", {
  set.seed(1)
  A <- matrix(stats::rnorm(16), 4); H <- crossprod(A) + diag(4)
  g <- stats::rnorm(4); theta <- stats::rnorm(4)
  pv <- list(lam = c(0, 2, 0.5, 3), alp = c(1, 1, 0.5, 1))
  b <- margEVT:::.cd_quadratic(theta, g, H, pv, lower = rep(-Inf, 4))
  obj <- function(x) sum(g * (x - theta)) +
    0.5 * drop(t(x - theta) %*% H %*% (x - theta)) +
    margEVT:::.penalty_exact(x, pv)
  opt <- stats::optim(theta, obj, method = "Nelder-Mead",
                      control = list(reltol = 1e-14, maxit = 20000))
  expect_lte(obj(b), opt$value + 1e-8)
  # a box bound is respected
  b2 <- margEVT:::.cd_quadratic(theta, g, H, pv, lower = c(-Inf, -Inf, -Inf, 5))
  expect_gte(b2[4], 5)
})

test_that("Newton and BFGS solvers give the same penalized fit", {
  df <- make_pen_df()
  args <- list(df, threshold = 16, loc_vars = c("x1", "x2", "x3"),
               scale_vars = c("x1", "x2"), penalty = "lasso", lambda = "bic",
               penalty_weights = "adaptive", obs_per_year = 50,
               verbose = FALSE)
  fb <- do.call(fit_nhpp, args)
  fn <- do.call(fit_nhpp, c(args, solver = "newton"))
  expect_true(fn$converged)
  expect_identical(fn$solver, "newton")
  expect_true(any(grepl("proximal Newton", capture.output(print(fn)))))
  expect_equal(fn$lambda, fb$lambda)
  expect_identical(fn$par_oper != 0, fb$par_oper != 0)
  expect_equal(fn$par_oper, fb$par_oper, tolerance = 1e-3)
  # exact zeros for the dropped coefficients
  expect_true(all(fn$par[fn$par_oper == 0] == 0))
  # the setting is inherited by refits
  expect_identical(refit_nhpp(fn, df)$solver, "newton")
})

test_that("the Newton solver gives up cleanly instead of looping", {
  df <- make_pen_df()
  fit <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2"),
                  penalty = "none", obs_per_year = 50, verbose = FALSE)
  dm <- fit$dm
  pv <- margEVT:::.penalty_vectors(dm, 0, 1, TRUE)
  bmap <- margEVT:::.shape_bound_map(dm, -0.5)
  pn <- function(init, ...) margEVT:::.prox_newton(init, dm, df$y, 16, pv, 50,
                                                  -0.5, bmap, ...)
  # infeasible start (exceedances beyond the upper endpoint)
  bad <- fit$par
  bad[length(bad)] <- -0.49
  bad[1] <- 0
  expect_null(pn(bad, tol_abs = 1e-6, active_tol = 1e-2))
  # KKT never satisfied: stops after max_iter
  expect_null(pn(fit$par + 0.1, tol_abs = -1, active_tol = 1e-2, max_iter = 2L))
  # a step that makes no progress
  testthat::local_mocked_bindings(.cd_quadratic = function(theta, ...) theta)
  expect_null(pn(fit$par + 0.1, tol_abs = -1, active_tol = 1e-2))
})

test_that("the Newton line search stops when no step decreases the objective", {
  df <- make_pen_df()
  fit <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2"),
                  penalty = "none", obs_per_year = 50, verbose = FALSE)
  dm <- fit$dm
  pv <- margEVT:::.penalty_vectors(dm, 0, 1, TRUE)
  calls <- 0L
  real <- margEVT:::.nhpp_eval
  testthat::local_mocked_bindings(.nhpp_eval = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) real(...) else list(value = 1e9, grad = NULL)
  })
  expect_null(margEVT:::.prox_newton(fit$par + 0.1, dm, df$y, 16, pv, 50,
                                     -Inf, NULL, -1, 1e-2))
  expect_gt(calls, 30L)
})

test_that("Newton handles the shape bound and falls back where needed", {
  set.seed(7L)
  n <- 3000L
  df <- data.frame(x = stats::rnorm(n))
  # bounded tail steeper than xi = -0.5: the bound is active
  df$y <- 10 + df$x + 3 * (1 - stats::runif(n)^(1 / 4))
  u <- unname(stats::quantile(df$y, 0.9))
  fb <- fit_nhpp(df, threshold = u, loc_vars = "x", penalty = "none",
                 obs_per_year = 100, verbose = FALSE)
  fn <- fit_nhpp(df, threshold = u, loc_vars = "x", penalty = "none",
                 obs_per_year = 100, solver = "newton", verbose = FALSE)
  expect_true(fn$converged)
  expect_equal(fn$shape_at_bound, fb$shape_at_bound)
  expect_equal(fn$par, fb$par, tolerance = 1e-3)
  # shape covariate under a finite bound: Newton does not apply, BFGS is used
  fs <- fit_nhpp(df, threshold = u, loc_vars = "x", shape_vars = "x",
                 penalty = "none", obs_per_year = 100, solver = "newton",
                 verbose = FALSE)
  expect_true(fs$converged)
})
