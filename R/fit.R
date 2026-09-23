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

.optim_nhpp <- function(...) {
  stats::optim(...)
}

# -----------------------------------------------------------------------------
# S3 constructor - internal
# -----------------------------------------------------------------------------
new_nhpp_fit <- function(par, par_oper, dm, threshold, nllh_pen, nllh_raw,
                         nllh_oper,
                         lambda, alpha, penalty, penalize_shape,
                         lambda_scaling, hessian, fitted, fitted_oper,
                         obs_per_year, active_tol,
                         converged, n_exc) {
  structure(
    list(
      par            = par,
      par_oper       = par_oper,
      dm             = dm,
      threshold      = threshold,
      nllh_pen       = nllh_pen,
      nllh_raw       = nllh_raw,
      nllh_oper      = nllh_oper,
      lambda         = lambda,
      alpha          = alpha,
      penalty        = penalty,
      penalize_shape = penalize_shape,
      lambda_scaling = lambda_scaling,
      hessian        = hessian,
      fitted         = fitted,
      fitted_oper    = fitted_oper,
      obs_per_year   = obs_per_year,
      active_tol     = active_tol,
      converged      = converged,
      n_exc          = n_exc
    ),
    class = "nhpp_fit"
  )
}

# -----------------------------------------------------------------------------
# .warmstart() - coordinate descent for a better starting point
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
# .fit_at_lambda() - internal workhorse
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
    .optim_nhpp(init, obj_fn, gr = gr_fn, method = "BFGS",
                control = list(maxit = maxit), hessian = calc_hessian),
    error = function(e) NULL
  )
  if (is.null(res) || res$convergence != 0)
    res <- tryCatch(
      .optim_nhpp(init, obj_fn, gr = gr_fn, method = "L-BFGS-B",
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
# fit_nhpp() - main user-facing function
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
#' @param active_tol Positive numeric tolerance used for BIC complexity
#'   counting, active-covariate extraction, and default coefficient reporting.
#'   Default \code{1e-2}.
#' @param lambda_scaling Character. Scaling used when \code{lambda = "bic"}:
#'   \code{"gradient"} calibrates the scale and shape penalties relative to
#'   the location-block score magnitude; \code{"common"} applies the same
#'   scalar lambda to all parameter blocks. Ignored for numeric \code{lambda}.
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
                     active_tol     = 1e-2,
                     lambda_scaling = c("gradient", "common"),
                     maxit          = 10000L,
                     calc_hessian   = FALSE,
                     verbose        = TRUE) {

  if (!is.data.frame(df))
    stop("fit_nhpp: `df` must be a data frame.")
  if (!"y" %in% names(df))
    stop("fit_nhpp: `df` must contain a column named `y` (the response).")
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold))
    stop("fit_nhpp: `threshold` must be a single numeric value that is finite.")
  if (!is.numeric(active_tol) || length(active_tol) != 1L ||
      !is.finite(active_tol) || active_tol <= 0)
    stop("fit_nhpp: `active_tol` must be a single positive numeric value.")
  if (!is.numeric(obs_per_year) || length(obs_per_year) != 1L ||
      !is.finite(obs_per_year) || obs_per_year <= 0)
    stop("fit_nhpp: `obs_per_year` must be a single positive numeric value.")
  if (!is.numeric(maxit) || length(maxit) != 1L ||
      !is.finite(maxit) || maxit < 1)
    stop("fit_nhpp: `maxit` must be a positive integer.")
  if (!is.logical(penalize_shape) || length(penalize_shape) != 1L ||
      is.na(penalize_shape))
    stop("fit_nhpp: `penalize_shape` must be TRUE or FALSE.")
  if (!is.logical(calc_hessian) || length(calc_hessian) != 1L ||
      is.na(calc_hessian))
    stop("fit_nhpp: `calc_hessian` must be TRUE or FALSE.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("fit_nhpp: `verbose` must be TRUE or FALSE.")

  penalty <- match.arg(penalty)
  lambda_scaling <- match.arg(lambda_scaling)

  y     <- df$y
  if (!is.numeric(y) || any(!is.finite(y)))
    stop("fit_nhpp: `df$y` must contain only finite numeric values.")
  n_exc <- sum(y > threshold, na.rm = TRUE)
  if (n_exc == 0L)
    stop("fit_nhpp: no observations exceed `threshold`.")
  if (n_exc < 5L && verbose)
    warning("fit_nhpp: fewer than 5 exceedances - estimates may be unreliable.")

  alpha <- switch(penalty,
                  none  = 0.5,
                  lasso = 1,
                  ridge = 0,
                  elnet = alpha
  )

  validate_block_control <- function(x, name, lower = -Inf, upper = Inf) {
    if (!is.numeric(x) || any(!is.finite(x)) ||
        length(x) < 1L || !length(x) %in% c(1L, 3L) ||
        any(x < lower) || any(x > upper))
      stop(sprintf(
        "fit_nhpp: `%s` must be a finite scalar or named vector c(mu=, sigma=, xi=).",
        name
      ))
    if (length(x) == 3L) {
      if (is.null(names(x)) || !setequal(names(x), c("mu", "sigma", "xi")))
        stop(sprintf(
          "fit_nhpp: a three-element `%s` vector must be named mu, sigma, and xi.",
          name
        ))
      x <- x[c("mu", "sigma", "xi")]
    }
    x
  }
  alpha <- validate_block_control(alpha, "alpha", lower = 0, upper = 1)

  if (penalty == "none") {
    lambda_resolved <- 0
  } else if (is.numeric(lambda)) {
    lambda_resolved <- validate_block_control(lambda, "lambda", lower = 0)
  } else if (identical(lambda, "bic")) {
    lambda_resolved <- NULL   # will be filled by grid search below
  } else {
    stop("fit_nhpp: `lambda` must be a non-negative numeric value or \"bic\".")
  }

  dm    <- build_design_matrices(df, loc_vars, scale_vars, shape_vars, free_vars)
  design_values <- unlist(dm[c("X_mu", "X_sigma", "X_xi")], use.names = FALSE)
  if (any(!is.finite(design_values)))
    stop("fit_nhpp: model covariates must contain only finite numeric values.")
  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

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
      active_tol     = active_tol,
      lambda_scaling = lambda_scaling,
      maxit          = maxit,
      verbose        = verbose
    )
    if (verbose)
      message(sprintf("fit_nhpp: selected lambda = %.5f", mean(lambda_resolved)))
  }

  if (any(lambda_resolved > 0))
    init <- .warmstart(init, dm, y, threshold, lambda_resolved, alpha,
                       penalize_shape, obs_per_year = obs_per_year)

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

  pen_idx_all <- c(
    dm$idx_pen_mu,
    p_mu + dm$idx_pen_sigma,
    if (penalize_shape) p_mu + p_sig + dm$idx_pen_xi else integer(0L)
  )
  par_oper <- par_hat
  small <- pen_idx_all[abs(par_oper[pen_idx_all]) < active_tol]
  if (length(small) > 0L) par_oper[small] <- 0
  nllh_oper <- pp_nllh(par_oper, dm, y, threshold,
                       lambda = 0, alpha = alpha,
                       pen_xi = penalize_shape,
                       obs_per_year = obs_per_year)

  beta_mu    <- par_hat[seq_len(p_mu)]
  beta_sigma <- par_hat[p_mu + seq_len(p_sig)]
  beta_xi    <- par_hat[p_mu + p_sig + seq_len(p_xi)]

  fitted <- list(
    mu    = as.numeric(dm$X_mu    %*% beta_mu),
    sigma = exp(as.numeric(dm$X_sigma %*% beta_sigma)),
    xi    = as.numeric(dm$X_xi    %*% beta_xi)
  )
  beta_mu_oper    <- par_oper[seq_len(p_mu)]
  beta_sigma_oper <- par_oper[p_mu + seq_len(p_sig)]
  beta_xi_oper    <- par_oper[p_mu + p_sig + seq_len(p_xi)]
  fitted_oper <- list(
    mu    = as.numeric(dm$X_mu    %*% beta_mu_oper),
    sigma = exp(as.numeric(dm$X_sigma %*% beta_sigma_oper)),
    xi    = as.numeric(dm$X_xi    %*% beta_xi_oper)
  )

  new_nhpp_fit(
    par            = par_hat,
    par_oper       = par_oper,
    dm             = dm,
    threshold      = threshold,
    nllh_pen       = res$nllh_pen,
    nllh_raw       = res$nllh_raw,
    nllh_oper      = nllh_oper,
    lambda         = lambda_resolved,
    alpha          = alpha,
    penalty        = penalty,
    penalize_shape = penalize_shape,
    lambda_scaling = lambda_scaling,
    hessian        = res$hessian,
    fitted         = fitted,
    fitted_oper    = fitted_oper,
    obs_per_year   = obs_per_year,
    active_tol     = active_tol,
    converged      = res$converged,
    n_exc          = n_exc
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.nhpp_fit <- function(x, ...) {
  cat("-- nhpp_fit --------------------------------------\n")
  cat(sprintf("  Threshold    : %.4g\n",  x$threshold))
  cat(sprintf("  Penalty      : %s\n",    x$penalty))
  cat(sprintf("  Lambda       : %.5g\n",  mean(x$lambda)))
  cat(sprintf("  Alpha        : %.3g\n",  mean(x$alpha)))
  cat(sprintf("  Converged    : %s\n",    x$converged))
  cat(sprintf("  nllh (raw)   : %.4f\n",  x$nllh_raw))
  cat(sprintf("  nllh (pen)   : %.4f\n",  x$nllh_pen))
  cat(sprintf("  obs/year     : %.2f\n",  x$obs_per_year))
  tol    <- if (is.null(x$active_tol)) 1e-2 else x$active_tol
  active <- x$par[.active_parameter_mask(x, tol)]
  cat(sprintf("  Active params: %d of %d\n", length(active), length(x$par)))
  cat("  Coefficients (non-zero):\n")
  print(round(active, 5L))
  invisible(x)
}

#' @export
coef.nhpp_fit <- function(object, operational = FALSE, ...) {
  if (isTRUE(operational) && !is.null(object$par_oper)) object$par_oper else
    object$par
}
