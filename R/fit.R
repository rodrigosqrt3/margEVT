# =============================================================================
# fit.R
# Core optimizer and S3 class constructor for nhpp_fit objects.
#
# Contents:
#   - new_nhpp_fit()     : S3 constructor (internal)
#   - .warmstart()       : coordinate descent warm start (internal)
#   - .fit_at_lambda()   : fit at a fixed lambda (internal workhorse)
#   - fit_nhpp()         : main user-facing fitting function
#   - print.nhpp_fit()   : print method
#   - coef.nhpp_fit()    : coef method
# =============================================================================

# -----------------------------------------------------------------------------
# S3 constructor — internal
# -----------------------------------------------------------------------------
new_nhpp_fit <- function(par, dm, threshold, nllh_pen, nllh_raw,
                         lambda, alpha, penalty, penalize_shape,
                         hessian, fitted, obs_per_year, converged) {
  structure(
    list(
      par            = par,
      dm             = dm,
      threshold      = threshold,
      nllh_pen       = nllh_pen,
      nllh_raw       = nllh_raw,
      lambda         = lambda,
      alpha          = alpha,
      penalty        = penalty,
      penalize_shape = penalize_shape,
      hessian        = hessian,
      fitted         = fitted,
      obs_per_year   = obs_per_year,
      converged      = converged
    ),
    class = "nhpp_fit"
  )
}

# -----------------------------------------------------------------------------
# .warmstart() — coordinate descent for a better starting point
# -----------------------------------------------------------------------------
.warmstart <- function(init, dm, y, threshold, lambda, alpha,
                       penalize_shape, n_passes = 5L, step = 0.01,
                       obs_per_year = 365.25) {

  if (all(lambda == 0)) return(init)

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  pen_idx_all <- c(
    dm$idx_pen_mu,
    p_mu + dm$idx_pen_sigma,
    if (penalize_shape) p_mu + p_sig + dm$idx_pen_xi else integer(0L)
  )
  if (length(pen_idx_all) == 0L) return(init)

  par     <- init
  obj     <- function(p) pp_nllh(p, dm, y, threshold, lambda, alpha,
                                 penalize_shape, obs_per_year = obs_per_year)
  f0      <- obj(par)
  lam_eff <- if (length(lambda) > 1L) mean(lambda) else lambda
  alp_eff <- if (length(alpha)  > 1L) mean(alpha)  else alpha

  for (pass in seq_len(n_passes)) {
    full_grad <- pp_grad(par, dm, y, threshold, lambda, alpha,
                         penalize_shape, obs_per_year = obs_per_year)
    for (j in pen_idx_all) {
      par_new    <- par
      par_new[j] <- par[j] - step * full_grad[j]
      par_new[j] <- sign(par_new[j]) *
        max(0, abs(par_new[j]) - lam_eff * alp_eff * step)
      f_new <- obj(par_new)
      if (is.finite(f_new) && f_new < f0) {
        par <- par_new
        f0  <- f_new
      }
    }
  }
  par
}

# -----------------------------------------------------------------------------
# .fit_at_lambda() — internal workhorse
# Fits the model at a fixed lambda. Called by fit_nhpp() and
# .select_lambda_bic(). Never called directly by the user.
# -----------------------------------------------------------------------------
.fit_at_lambda <- function(dm, y, threshold, lambda, alpha,
                           penalize_shape, init = NULL,
                           maxit = 10000L, calc_hessian = FALSE,
                           obs_per_year = 365.25) {

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  if (is.null(init)) {
    y_exc   <- y[y > attr(dm, "threshold")]
    init    <- rep(0, p_mu + p_sig + p_xi)
    if (length(y_exc) > 0L) {
      init[1L]        <- mean(y_exc)
      init[p_mu + 1L] <- log(max(stats::sd(y_exc), 0.1))
    }
  }

  obj_fn <- function(par) pp_nllh(par, dm, y, threshold, lambda, alpha,
                                  penalize_shape, obs_per_year = obs_per_year)
  gr_fn  <- function(par) pp_grad(par, dm, y, threshold, lambda, alpha,
                                  penalize_shape, obs_per_year = obs_per_year)

  res <- tryCatch(
    stats::optim(init, obj_fn, gr = gr_fn, method = "BFGS",
                 control = list(maxit = maxit), hessian = calc_hessian),
    error = function(e) NULL
  )
  if (is.null(res) || res$convergence != 0)
    res <- tryCatch(
      stats::optim(init, obj_fn, gr = gr_fn, method = "L-BFGS-B",
                   control = list(maxit = maxit, factr = 1e7),
                   hessian = calc_hessian),
      error = function(e) NULL
    )

  if (is.null(res))
    return(list(converged = FALSE, par = init,
                nllh_pen = NA_real_, nllh_raw = NA_real_,
                hessian = NULL))

  par_hat      <- res$par
  nllh_raw     <- pp_nllh(par_hat, dm, y, threshold,
                          lambda = 0, alpha = alpha,
                          pen_xi = penalize_shape,
                          obs_per_year = obs_per_year)

  list(
    converged = (res$convergence == 0L),
    par       = par_hat,
    nllh_pen  = res$value,
    nllh_raw  = nllh_raw,
    hessian   = if (calc_hessian) res$hessian else NULL
  )
}

# -----------------------------------------------------------------------------
# fit_nhpp() — main user-facing function
# -----------------------------------------------------------------------------

#' Fit a non-homogeneous point process model for extremes
#'
#' Fits a peaks-over-threshold point process model via penalized maximum
#' likelihood with BFGS optimization and exact analytical gradients.
#'
#' @param df A data frame containing the response and all covariates.
#'   Must have a column named \code{y} (the observations).
#' @param threshold Numeric scalar. The extreme value threshold u.
#' @param loc_vars Character vector of covariate names for the location
#'   parameter. \code{NULL} for stationary location.
#' @param scale_vars Character vector of covariate names for the scale
#'   parameter. \code{NULL} for stationary scale.
#' @param shape_vars Character vector of covariate names for the shape
#'   parameter. \code{NULL} for stationary shape.
#' @param free_vars Character vector of covariate names that are never
#'   penalized (e.g. seasonality terms). \code{NULL} by default.
#' @param penalty Character. Penalty type: \code{"none"} (pure MLE),
#'   \code{"lasso"} (L1), \code{"ridge"} (L2), or \code{"elnet"}
#'   (elastic net, requires \code{alpha}).
#' @param alpha Numeric in [0, 1]. Elastic-net mixing parameter.
#'   Only used when \code{penalty = "elnet"}. 1 = LASSO, 0 = ridge.
#'   Scalar or named vector \code{c(mu = , sigma = , xi = )}.
#' @param lambda Numeric scalar or named vector \code{c(mu=, sigma=, xi=)},
#'   or \code{"bic"} to select automatically via a two-phase BIC grid search.
#'   Ignored when \code{penalty = "none"}.
#' @param penalize_shape Logical. Penalize shape parameter covariates?
#'   Default \code{TRUE}.
#' @param obs_per_year Numeric. Observations per year. E.g. \code{365.25}
#'   for daily data, \code{52} for weekly. Default \code{365.25}.
#' @param maxit Integer. Maximum optimizer iterations. Default \code{10000L}.
#' @param calc_hessian Logical. Compute Hessian at solution? Needed for
#'   delta-method standard errors. Default \code{FALSE}.
#' @param verbose Logical. Print progress during BIC grid search and
#'   convergence warnings. Default \code{TRUE}.
#'
#' @return An object of class \code{nhpp_fit}.
#'
#' @export
fit_nhpp <- function(df, threshold,
                     loc_vars       = NULL,
                     scale_vars     = NULL,
                     shape_vars     = NULL,
                     free_vars      = NULL,
                     penalty        = c("none", "lasso", "ridge", "elnet"),
                     alpha          = 0.5,
                     lambda         = "bic",
                     penalize_shape = TRUE,
                     obs_per_year   = 365.25,
                     maxit          = 10000L,
                     calc_hessian   = FALSE,
                     verbose        = TRUE) {

  # ── Validate inputs ────────────────────────────────────────────────────────
  if (!is.data.frame(df))
    stop("fit_nhpp: `df` must be a data frame.")
  if (!"y" %in% names(df))
    stop("fit_nhpp: `df` must contain a column named `y` (the response).")
  if (!is.numeric(threshold) || length(threshold) != 1L)
    stop("fit_nhpp: `threshold` must be a single numeric value.")

  penalty <- match.arg(penalty)

  y     <- df$y
  n_exc <- sum(y > threshold, na.rm = TRUE)
  if (n_exc < 5L && verbose)
    warning("fit_nhpp: fewer than 5 exceedances — estimates may be unreliable.")

  # ── Resolve alpha from penalty type ───────────────────────────────────────
  alpha <- switch(penalty,
                  none  = 0.5,   # irrelevant, lambda will be 0
                  lasso = 1,
                  ridge = 0,
                  elnet = alpha
  )

  # ── Resolve lambda from penalty type ──────────────────────────────────────
  if (penalty == "none") {
    lambda_resolved <- 0
  } else if (is.numeric(lambda)) {
    lambda_resolved <- lambda
  } else if (identical(lambda, "bic")) {
    lambda_resolved <- NULL   # will be filled by grid search below
  } else {
    stop("fit_nhpp: `lambda` must be a positive numeric value or \"bic\".")
  }

  # ── Build design matrices ──────────────────────────────────────────────────
  dm    <- build_design_matrices(df, loc_vars, scale_vars, shape_vars, free_vars)
  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  # ── Initial parameter vector ───────────────────────────────────────────────
  y_exc <- y[y > threshold]
  init  <- rep(0, p_mu + p_sig + p_xi)
  if (length(y_exc) > 0L) {
    init[1L]        <- mean(y_exc)
    init[p_mu + 1L] <- log(max(stats::sd(y_exc), 0.1))
  }
  names(init) <- c(
    paste0("mu.",    colnames(dm$X_mu)),
    paste0("sigma.", colnames(dm$X_sigma)),
    paste0("xi.",    colnames(dm$X_xi))
  )

  # ── BIC grid search if needed ──────────────────────────────────────────────
  if (is.null(lambda_resolved)) {
    if (verbose)
      message("fit_nhpp: running BIC lambda selection (penalty = '",
              penalty, "')...")
    lambda_resolved <- .select_lambda_bic(
      dm             = dm,
      y              = y,
      threshold      = threshold,
      alpha          = alpha,
      penalize_shape = penalize_shape,
      init           = init,
      obs_per_year   = obs_per_year,
      maxit          = maxit,
      verbose        = verbose
    )
    if (verbose)
      message(sprintf("fit_nhpp: selected lambda = %.5f", mean(lambda_resolved)))
  }

  # ── Warm start ────────────────────────────────────────────────────────────
  if (any(lambda_resolved > 0))
    init <- .warmstart(init, dm, y, threshold, lambda_resolved, alpha,
                       penalize_shape, obs_per_year = obs_per_year)

  # ── Fit at resolved lambda ─────────────────────────────────────────────────
  res <- .fit_at_lambda(dm, y, threshold,
                        lambda       = lambda_resolved,
                        alpha        = alpha,
                        penalize_shape = penalize_shape,
                        init         = init,
                        maxit        = maxit,
                        calc_hessian = calc_hessian,
                        obs_per_year = obs_per_year)

  if (!res$converged && verbose)
    warning("fit_nhpp: optimizer did not converge.")

  par_hat <- res$par
  names(par_hat) <- names(init)

  # ── Fitted values ──────────────────────────────────────────────────────────
  beta_mu    <- par_hat[seq_len(p_mu)]
  beta_sigma <- par_hat[p_mu + seq_len(p_sig)]
  beta_xi    <- par_hat[p_mu + p_sig + seq_len(p_xi)]

  fitted <- list(
    mu    = as.numeric(dm$X_mu    %*% beta_mu),
    sigma = exp(as.numeric(dm$X_sigma %*% beta_sigma)),
    xi    = as.numeric(dm$X_xi    %*% beta_xi)
  )

  new_nhpp_fit(
    par            = par_hat,
    dm             = dm,
    threshold      = threshold,
    nllh_pen       = res$nllh_pen,
    nllh_raw       = res$nllh_raw,
    lambda         = lambda_resolved,
    alpha          = alpha,
    penalty        = penalty,
    penalize_shape = penalize_shape,
    hessian        = res$hessian,
    fitted         = fitted,
    obs_per_year   = obs_per_year,
    converged      = res$converged
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.nhpp_fit <- function(x, ...) {
  cat("── nhpp_fit ──────────────────────────────────────\n")
  cat(sprintf("  Threshold    : %.4g\n",  x$threshold))
  cat(sprintf("  Penalty      : %s\n",    x$penalty))
  cat(sprintf("  Lambda       : %.5g\n",  mean(x$lambda)))
  cat(sprintf("  Alpha        : %.3g\n",  mean(x$alpha)))
  cat(sprintf("  Converged    : %s\n",    x$converged))
  cat(sprintf("  nllh (raw)   : %.4f\n",  x$nllh_raw))
  cat(sprintf("  nllh (pen)   : %.4f\n",  x$nllh_pen))
  cat(sprintf("  obs/year     : %.2f\n",  x$obs_per_year))
  active <- x$par[abs(x$par) > 1e-4]
  cat(sprintf("  Active params: %d of %d\n", length(active), length(x$par)))
  cat("  Coefficients (non-zero):\n")
  print(round(active, 5L))
  invisible(x)
}

#' @export
coef.nhpp_fit <- function(object, ...) {
  object$par
}
