# -----------------------------------------------------------------------------
# Test Helpers
# -----------------------------------------------------------------------------
make_select_df <- function(n = 400L, seed = 7L) {
  set.seed(seed)
  data.frame(
    y = c(stats::rexp(n, 0.3), stats::runif(40L, 5, 20)),
    x = stats::rnorm(n + 40L),
    z = stats::rnorm(n + 40L)
  )
}

# -----------------------------------------------------------------------------
# Core Functionality Tests
# -----------------------------------------------------------------------------
test_that("lambda=bic returns a converged fit", {
  df  <- make_select_df()
  expect_message(
    fit <- fit_nhpp(df, threshold = 4,
                    loc_vars = c("x", "z"),
                    penalty  = "lasso",
                    lambda   = "bic",
                    verbose  = TRUE),
    # The BIC is flat over a range of lambdas here, so which branch is taken
    # depends on platform-specific rounding.
    regexp = "Fine search improved|Coarse optimum confirmed"
  )
  expect_s3_class(fit, "nhpp_fit")
  expect_true(fit$converged)
  expect_true(is.numeric(fit$lambda))
  expect_true(all(fit$lambda > 0))
  expect_equal(fit$active_tol, 1e-2)
})

test_that("bic selection returns a named lambda vector", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = "x",
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_named(fit$lambda, c("mu", "sigma", "xi"))
})

test_that("bic lasso shrinks noise covariate more than signal", {
  set.seed(99L)
  n   <- 500L
  y   <- c(stats::rexp(n, 0.3), stats::runif(50L, 5, 20))
  x   <- c(stats::rnorm(n), stats::rnorm(50L, mean = 2))
  z   <- stats::rnorm(n + 50L)
  df  <- data.frame(y = y, x = x, z = z)
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = c("x", "z"),
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_gte(abs(fit$par["mu.x"]), abs(fit$par["mu.z"]))
})

test_that("ridge bic fit has non-zero coefficients for all covariates", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = c("x", "z"),
                  penalty  = "ridge",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_gt(abs(fit$par["mu.x"]), 1e-6)
  expect_gt(abs(fit$par["mu.z"]), 1e-6)
})

# -----------------------------------------------------------------------------
# Coverage & Branch Tests
# -----------------------------------------------------------------------------
test_that("bic selection handles no location covariates to penalize", {
  df  <- make_select_df()
  fit <- fit_nhpp(df, threshold = 4,
                  loc_vars = NULL,
                  scale_vars = "x",
                  penalty  = "lasso",
                  lambda   = "bic",
                  verbose  = FALSE)
  expect_s3_class(fit, "nhpp_fit")
})

test_that("fit_nhpp rejects non-finite model covariates before BIC selection", {
  df  <- make_select_df()
  df$x[1] <- NaN
  expect_error(
    fit_nhpp(df, threshold = 4,
             loc_vars = "x",
             penalty  = "lasso",
             lambda   = "bic",
             verbose  = FALSE),
    regexp = "model covariates must contain only finite numeric values"
  )
})

test_that("bic selection throws error when all coarse grid fits fail", {
  df <- make_select_df()
  dm <- margEVT:::build_design_matrices(
    df, loc_vars = "x", scale_vars = NULL,
    shape_vars = NULL, free_vars = NULL
  )
  init <- rep(0, ncol(dm$X_mu) + ncol(dm$X_sigma) + ncol(dm$X_xi))

  testthat::local_mocked_bindings(
    .fit_at_lambda = function(...) list(
      converged = FALSE,
      par = init,
      nllh_pen = NA_real_,
      nllh_raw = NA_real_,
      hessian = NULL
    )
  )

  expect_error(
    margEVT:::.select_lambda_bic(
      dm = dm, y = df$y, threshold = 4, alpha = 1,
      penalize_shape = TRUE, init = init, verbose = FALSE
    ),
    regexp = "all coarse grid fits failed"
  )
})

test_that("bic selection rejects fits with invalid operational likelihood", {
  df <- make_select_df()
  dm <- margEVT:::build_design_matrices(
    df, loc_vars = "x", scale_vars = NULL,
    shape_vars = NULL, free_vars = NULL
  )
  init <- rep(0, ncol(dm$X_mu) + ncol(dm$X_sigma) + ncol(dm$X_xi))

  testthat::local_mocked_bindings(
    .fit_at_lambda = function(...) list(
      converged = TRUE,
      par = init,
      nllh_pen = 0,
      nllh_raw = 0,
      hessian = NULL
    ),
    pp_nllh = function(...) Inf
  )

  expect_error(
    margEVT:::.select_lambda_bic(
      dm = dm, y = df$y, threshold = 4, alpha = 1,
      penalize_shape = TRUE, init = init, verbose = FALSE
    ),
    regexp = "all coarse grid fits failed"
  )
})

test_that("bic selection stops early when valley is found (Mocked via testthat 3)", {
  s <- make_select_df()
  dm <- margEVT:::build_design_matrices(s, loc_vars = "x", scale_vars = NULL, shape_vars = NULL, free_vars = NULL)
  init <- rep(0, ncol(dm$X_mu) + ncol(dm$X_sigma) + ncol(dm$X_xi))

  counter <- 0
  mock_fit <- function(...) {
    counter <<- counter + 1

    list(
      converged = TRUE,
      # Encode the path position in the fitted vector. Since selection now
      # recomputes the unpenalized likelihood at the thresholded operational
      # estimator, the likelihood mock below reads this value.
      par = rep(counter, length(init)),
      nllh_pen = 0,
      nllh_raw = 0,
      hessian = NULL
    )
  }

  mock_nllh <- function(par, ...) 100 + unname(par[1]) * 40

  testthat::local_mocked_bindings(
    .fit_at_lambda = mock_fit,
    pp_nllh = mock_nllh
  )

  expect_message(
    fit <- margEVT:::.select_lambda_bic(
      dm = dm, y = s$y, threshold = 4, alpha = 1,
      penalize_shape = TRUE, init = init, verbose = TRUE
    ),
    regexp = "Valley found, stopping coarse search early"
  )
})

test_that("the fine grid improves on the coarse optimum (Mocked)", {
  s <- make_select_df()
  dm <- margEVT:::build_design_matrices(s, loc_vars = "x", scale_vars = NULL,
                                        shape_vars = NULL, free_vars = NULL)
  init <- rep(0, ncol(dm$X_mu) + ncol(dm$X_sigma) + ncol(dm$X_xi))

  # The fitted vector encodes log(lambda); the likelihood is minimized at
  # lambda = 3.3, which lies strictly between two coarse grid points.
  testthat::local_mocked_bindings(
    .fit_at_lambda = function(..., lambda) list(
      converged = TRUE,
      par = rep(log(lambda[1L]), length(init)),
      nllh_pen = 0, nllh_raw = 0, hessian = NULL
    ),
    pp_nllh = function(par, ...) 1000 * (unname(par[1L]) - log(3.3))^2
  )

  expect_message(
    sel <- margEVT:::.select_lambda_bic(
      dm = dm, y = s$y, threshold = 4, alpha = 1,
      penalize_shape = TRUE, init = init, verbose = TRUE
    ),
    regexp = "Fine search improved"
  )
  # strictly between the neighbouring coarse points 2.944 and 3.864
  expect_gt(unname(sel$lambda["mu"]), 2.95)
  expect_lt(unname(sel$lambda["mu"]), 3.86)
})

test_that("gradient lambda scaling calibrates block penalties", {
  df <- make_select_df()
  grid <- c(20, 5, 1)
  fit_all <- fit_nhpp(df, threshold = 4, loc_vars = "x", scale_vars = "x",
                      shape_vars = "z", penalty = "lasso", lambda = "bic",
                      lambda_scaling = "gradient", lambda_grid = grid,
                      verbose = FALSE)
  expect_identical(fit_all$lambda_scaling, "gradient")
  expect_false(isTRUE(all.equal(unname(fit_all$lambda["mu"]),
                                unname(fit_all$lambda["sigma"]))))
  # No penalized location coefficients: the location score scale falls back
  # to one, and missing blocks inherit the location scale.
  fit_noloc <- fit_nhpp(df, threshold = 4, loc_vars = NULL, scale_vars = "x",
                        penalty = "lasso", lambda = "bic",
                        lambda_scaling = "gradient", lambda_grid = grid,
                        verbose = FALSE)
  expect_equal(unname(fit_noloc$lambda["xi"]), unname(fit_noloc$lambda["mu"]))
  fit_loc <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "lasso",
                      lambda = "bic", lambda_scaling = "gradient",
                      lambda_grid = grid, verbose = FALSE)
  expect_equal(unname(fit_loc$lambda), rep(unname(fit_loc$lambda["mu"]), 3L))
})
