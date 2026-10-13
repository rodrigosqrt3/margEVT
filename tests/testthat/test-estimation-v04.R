# Tests for estimation features introduced in margEVT 0.4.0:
# coordinate penalty vectors, KKT convergence, adaptive and numeric penalty
# weights, user starting values and grids, path storage, and the fact that
# the returned model is the BIC-scored path solution.

make_v04_df <- function(n = 1500L, seed = 4L) {
  set.seed(seed)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   x3 = stats::rnorm(n))
  df$y <- 10 + 2 * df$x1 + stats::rexp(n, 1 / 3)
  df
}

test_that(".resolve_block handles scalar, named and unnamed vectors", {
  expect_equal(margEVT:::.resolve_block(2, "sigma"), 2)
  expect_equal(margEVT:::.resolve_block(c(xi = 3, mu = 1, sigma = 2), "sigma"), 2)
  expect_equal(margEVT:::.resolve_block(c(1, 2, 3), "xi"), 3)
})

test_that(".penalty_vectors assigns block lambdas and weights by coordinate", {
  df <- make_v04_df(200L)
  dm <- build_design_matrices(df, loc_vars = "x1", scale_vars = "x2",
                              shape_vars = "x3")
  pv <- margEVT:::.penalty_vectors(dm, c(mu = 1, sigma = 2, xi = 3), 1, TRUE)
  expect_equal(pv$lam, c(0, 1, 0, 2, 0, 3))
  pv2 <- margEVT:::.penalty_vectors(dm, 1, 0.5, FALSE,
                                    penalty_weights = c(9, 2, 9, 3, 9, 4))
  expect_equal(pv2$lam, c(0, 2, 0, 3, 0, 0))
  expect_equal(pv2$alp, c(0, 0.5, 0, 0.5, 0, 0))
})

test_that("weighted pp_nllh and pp_grad agree with finite differences", {
  df <- make_v04_df(400L)
  dm <- build_design_matrices(df, loc_vars = c("x1", "x2"))
  par <- c(15, 1, 0.2, log(3), 0.05)
  w <- c(0, 2, 0.5, 0, 0)
  g <- pp_grad(par, dm, df$y, 16, lambda = 1.5, alpha = 0.7,
               penalty_weights = w)
  fd <- vapply(seq_along(par), function(k) {
    h <- 1e-6; a <- b <- par; a[k] <- a[k] + h; b[k] <- b[k] - h
    (pp_nllh(a, dm, df$y, 16, 1.5, 0.7, penalty_weights = w) -
       pp_nllh(b, dm, df$y, 16, 1.5, 0.7, penalty_weights = w)) / (2 * h)
  }, numeric(1L))
  expect_equal(g, fd, tolerance = 1e-4)
})

test_that(".objective_factory caches and flags the support barrier", {
  df <- make_v04_df(300L)
  dm <- build_design_matrices(df)
  pv <- margEVT:::.penalty_vectors(dm, 0, 1, TRUE)
  obj <- margEVT:::.objective_factory(dm, df$y, 16, pv)
  par <- c(15, log(3), 0.05)
  expect_identical(obj$fn(par), obj$fn(par))
  expect_length(obj$gr(par), 3L)
  expect_equal(obj$fn(c(15, log(3), -5)), 1e9)
})

test_that(".kkt_violation uses the subgradient at zero and rejects non-finite input", {
  pv <- list(lam = c(0, 2), alp = c(0, 1))
  expect_equal(margEVT:::.kkt_violation(c(1, 0), c(0.1, 1.5), pv, 1e-2), 0.1)
  expect_equal(margEVT:::.kkt_violation(c(1, 0), c(0, 3), pv, 1e-2), 1)
  expect_identical(margEVT:::.kkt_violation(c(1, 0), c(NA, 0), pv, 1e-2), Inf)
})

test_that(".warmstart applies coordinate-specific soft thresholds", {
  df <- make_v04_df(500L)
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

test_that(".fit_at_lambda tries L-BFGS-B and a BFGS restart when KKT fails", {
  df <- make_v04_df(300L)
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
  expect_identical(calls, c("BFGS", "L-BFGS-B", "BFGS"))
  expect_false(res$converged)
  expect_true(is.finite(res$kkt))
})

test_that(".fit_at_lambda discards optimizer results outside the support", {
  df <- make_v04_df(300L)
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

test_that("BIC fit returns the path solution scored by the criterion", {
  df <- make_v04_df()
  fit <- fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2", "x3"),
                  penalty = "lasso", lambda = "bic", verbose = FALSE)
  expect_true(fit$converged)
  tab <- fit$path$table
  expect_true(all(c("lambda", "phase", "bic", "k_active", "converged",
                    "selected") %in% names(tab)))
  expect_equal(sum(tab$selected), 1L)
  sel <- which(tab$selected)
  expect_equal(unname(fit$par_oper), unname(fit$path$coef[sel, ]),
               tolerance = 1e-5)
  expect_equal(bic_nhpp(fit), tab$bic[sel], tolerance = 1e-6)
  expect_equal(unname(fit$lambda["mu"]), tab$lambda[sel])
})

test_that("user lambda_grid is evaluated in one pass and validated", {
  df <- make_v04_df()
  fit <- suppressMessages(
    fit_nhpp(df, threshold = 16, loc_vars = c("x1", "x2"),
             penalty = "lasso", lambda = "bic",
             lambda_grid = c(1, 50, 5), verbose = TRUE))
  expect_identical(unique(fit$path$table$phase), "user")
  expect_equal(fit$path$table$lambda, c(50, 5, 1))
  expect_error(fit_nhpp(df, threshold = 16, loc_vars = "x1",
                        penalty = "lasso", lambda_grid = 1, verbose = FALSE),
               "lambda_grid")
  dm <- build_design_matrices(df, loc_vars = "x1")
  testthat::local_mocked_bindings(
    .fit_at_lambda = function(...) list(converged = FALSE, par = NULL,
                                        nllh_pen = NA, nllh_raw = NA)
  )
  expect_error(
    margEVT:::.select_lambda_bic(dm, df$y, 16, 1, TRUE, c(15, 0, 1, 0),
                                 lambda_grid = c(2, 1), verbose = FALSE),
    "all grid fits failed")
})

test_that("adaptive weights are stored and applied", {
  df <- make_v04_df()
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
  df <- make_v04_df(400L)
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
  df <- make_v04_df(600L)
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
  df <- make_v04_df(400L)
  fit <- fit_nhpp(df, threshold = 16, loc_vars = "x1", penalty = "none",
                  penalty_weights = "adaptive", verbose = FALSE)
  expect_null(fit$penalty_weights)
  expect_identical(fit$weights_source, "none")
})
