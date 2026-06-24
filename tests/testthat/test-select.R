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
    regexp = "Fine search improved"
  )
  expect_s3_class(fit, "nhpp_fit")
  expect_true(fit$converged)
  expect_true(is.numeric(fit$lambda))
  expect_true(all(fit$lambda > 0))
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

test_that("bic selection throws error when all coarse grid fits fail", {
  df  <- make_select_df()
  df$x[1] <- NaN
  expect_error(
    fit_nhpp(df, threshold = 4,
             loc_vars = "x",
             penalty  = "lasso",
             lambda   = "bic",
             verbose  = FALSE),
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
    nllh <- 100 + counter * 40

    list(
      converged = TRUE,
      par = rep(1, length(init)),
      nllh_pen = nllh,
      nllh_raw = nllh,
      hessian = NULL
    )
  }

  testthat::local_mocked_bindings(.fit_at_lambda = mock_fit)

  expect_message(
    fit <- margEVT:::.select_lambda_bic(
      dm = dm, y = s$y, threshold = 4, alpha = 1,
      penalize_shape = TRUE, init = init, verbose = TRUE
    ),
    regexp = "Valley found, stopping coarse search early"
  )
})
