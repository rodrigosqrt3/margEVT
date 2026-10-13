# =============================================================================
# fit.R
# Core optimizer and S3 class constructor for nhpp_fit objects.
#
# Contents:
#   - new_nhpp_fit()       : S3 constructor (internal)
#   - .penalized_index()   : penalized coordinates of a design (internal)
#   - .warmstart()         : proximal coordinate warm start (internal)
#   - .fit_at_lambda()     : fit at a fixed lambda (internal workhorse)
#   - .adaptive_weights()  : adaptive-LASSO weights (internal)
#   - fit_nhpp()           : main user-facing fitting function
#   - print.nhpp_fit(), coef.nhpp_fit()
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
                         converged, n_exc,
                         penalty_weights = NULL, weights_source = "none",
                         kkt = NA_real_, path = NULL, y = NULL,
                         grad_tol = 5e-4, adaptive_gamma = 1) {
  structure(
    list(
      par             = par,
      par_oper        = par_oper,
      dm              = dm,
      threshold       = threshold,
      nllh_pen        = nllh_pen,
      nllh_raw        = nllh_raw,
      nllh_oper       = nllh_oper,
      lambda          = lambda,
      alpha           = alpha,
      penalty         = penalty,
      penalize_shape  = penalize_shape,
      lambda_scaling  = lambda_scaling,
      hessian         = hessian,
      fitted          = fitted,
      fitted_oper     = fitted_oper,
      obs_per_year    = obs_per_year,
      active_tol      = active_tol,
      converged       = converged,
      n_exc           = n_exc,
      penalty_weights = penalty_weights,
      weights_source  = weights_source,
      kkt             = kkt,
      path            = path,
      y               = y,
      grad_tol        = grad_tol,
      adaptive_gamma  = adaptive_gamma
    ),
    class = "nhpp_fit"
  )
}

# Penalized coordinates of the full parameter vector.
.penalized_index <- function(dm, penalize_shape) {
  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  c(dm$idx_pen_mu,
    p_mu + dm$idx_pen_sigma,
    if (isTRUE(penalize_shape)) p_mu + p_sig + dm$idx_pen_xi else integer(0L))
}

# -----------------------------------------------------------------------------
# .warmstart() - proximal coordinate warm start
# Soft thresholds use the coordinate-specific penalty lambda_k * alpha_k, so
# block-specific and weighted penalties are respected.
# -----------------------------------------------------------------------------
.warmstart <- function(init, dm, y, threshold, lambda, alpha,
                       penalize_shape, n_passes = 5L, step = 0.01,
                       obs_per_year = 365.25, penalty_weights = NULL) {

  if (all(lambda == 0)) return(init)

  pv  <- .penalty_vectors(dm, lambda, alpha, penalize_shape, penalty_weights)
  idx <- which(pv$lam > 0)
  if (length(idx) == 0L) return(init)

  obj <- .objective_factory(dm, y, threshold, pv, obs_per_year = obs_per_year)
  par <- init
  f0  <- obj$fn(par)

  for (pass in seq_len(n_passes)) {
    full_grad <- obj$gr(par)
    for (j in idx) {
      par_new    <- par
      par_new[j] <- par[j] - step * full_grad[j]
      par_new[j] <- sign(par_new[j]) *
        max(0, abs(par_new[j]) - pv$lam[j] * pv$alp[j] * step)
      f_new <- obj$fn(par_new)
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
# Fits the model at a fixed lambda. BFGS is tried first; if the KKT check
# fails, L-BFGS-B and a final BFGS restart are attempted. Convergence is
# declared only when the KKT violation (scaled by the exceedance count) is
# below `grad_tol`.
# -----------------------------------------------------------------------------
.fit_at_lambda <- function(dm, y, threshold, lambda, alpha,
                           penalize_shape, init = NULL,
                           maxit = 10000L, calc_hessian = FALSE,
                           obs_per_year = 365.25, penalty_weights = NULL,
                           grad_tol = 5e-4, active_tol = 1e-2) {

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  if (is.null(init)) {
    y_exc <- y[y > threshold]
    init  <- rep(0, p_mu + p_sig + p_xi)
    if (length(y_exc) > 0L) {
      init[1L]        <- mean(y_exc)
      init[p_mu + 1L] <- log(max(stats::sd(y_exc), 0.1))
    }
  }

  pv  <- .penalty_vectors(dm, lambda, alpha, penalize_shape, penalty_weights)
  obj <- .objective_factory(dm, y, threshold, pv, obs_per_year = obs_per_year)
  n_exc   <- sum(y > threshold)
  tol_abs <- grad_tol * max(1, n_exc)

  run <- function(start, method) {
    control <- if (method == "BFGS")
      list(maxit = maxit, reltol = 1e-10) else
        list(maxit = maxit, factr = 1e7, pgtol = 0)
    res <- tryCatch(
      .optim_nhpp(start, obj$fn, gr = obj$gr, method = method,
                  control = control),
      error = function(e) NULL)
    if (is.null(res) || any(!is.finite(res$par)) || !is.finite(res$value) ||
        res$value >= 1e9) return(NULL)
    res$kkt <- .kkt_violation(res$par, obj$gr(res$par), pv, active_tol)
    res
  }

  ok <- function(res) !is.null(res) && res$kkt <= tol_abs
  res1 <- run(init, "BFGS")
  res2 <- res3 <- NULL
  if (!ok(res1)) {
    res2 <- run(if (!is.null(res1)) res1$par else init, "L-BFGS-B")
    if (!ok(res2) && !is.null(res2))
      res3 <- run(res2$par, "BFGS")
  }
  attempts <- Filter(Negate(is.null), list(res1, res2, res3))

  if (length(attempts) == 0L)
    return(list(converged = FALSE, par = init,
                nllh_pen = NA_real_, nllh_raw = NA_real_,
                hessian = NULL, kkt = NA_real_))

  kkts <- vapply(attempts, `[[`, numeric(1L), "kkt")
  best <- attempts[[which.min(kkts)]]
  par_hat <- best$par

  nllh_raw <- pp_nllh(par_hat, dm, y, threshold,
                      lambda = 0, alpha = alpha,
                      pen_xi = penalize_shape,
                      obs_per_year = obs_per_year)

  list(
    converged = best$kkt <= tol_abs,
    par       = par_hat,
    nllh_pen  = best$value,
    nllh_raw  = nllh_raw,
    hessian   = if (calc_hessian)
      stats::optimHess(par_hat, obj$fn, obj$gr) else NULL,
    kkt       = best$kkt / max(1, n_exc)
  )
}

# -----------------------------------------------------------------------------
# .adaptive_weights() - adaptive-LASSO weights w_k = 1 / |theta_init_k|^gamma
# The initial estimate is the unpenalized MLE; if it fails the KKT check, a
# light ridge fit (alpha = 0) is used with increasing strength.
# -----------------------------------------------------------------------------
.adaptive_weights <- function(dm, y, threshold, init, penalize_shape,
                              obs_per_year, gamma = 1, floor_abs = 1e-4,
                              maxit = 10000L, grad_tol = 5e-4,
                              active_tol = 1e-2,
                              ridge_lambdas = c(0.1, 1, 10)) {
  pen_idx <- .penalized_index(dm, penalize_shape)
  base <- .fit_at_lambda(dm, y, threshold, lambda = 0, alpha = 1,
                         penalize_shape = penalize_shape, init = init,
                         maxit = maxit, obs_per_year = obs_per_year,
                         grad_tol = grad_tol, active_tol = active_tol)
  source <- "mle"
  if (!isTRUE(base$converged)) {
    source <- NA_character_
    for (lr in ridge_lambdas) {
      base <- .fit_at_lambda(dm, y, threshold, lambda = lr, alpha = 0,
                             penalize_shape = penalize_shape, init = init,
                             maxit = maxit, obs_per_year = obs_per_year,
                             grad_tol = grad_tol, active_tol = active_tol)
      if (isTRUE(base$converged)) {
        source <- paste0("ridge(lambda=", lr, ")")
        break
      }
    }
  }
  if (is.na(source))
    stop("fit_nhpp: adaptive weights need a convergent initial fit; ",
         "supply numeric `penalty_weights` instead.")
  w <- numeric(length(init))
  w[pen_idx] <- 1 / pmax(abs(base$par[pen_idx]), floor_abs)^gamma
  names(w) <- names(init)
  list(weights = w, source = source)
}

# -----------------------------------------------------------------------------
# fit_nhpp() - main user-facing function
# -----------------------------------------------------------------------------

#' Fit a non-homogeneous point process model for extremes
#'
#' Fits a peaks-over-threshold point process model via penalized maximum
#' likelihood with quasi-Newton optimization and exact analytical gradients.
#' The location, log-scale, and shape parameters are linear in the supplied
#' covariates; seasonal harmonics and other \code{free_vars} are never
#' penalized.
#'
#' Convergence is assessed with the Karush-Kuhn-Tucker (KKT) conditions of
#' the smoothed elastic-net objective: unpenalized and active coordinates
#' must have a near-zero gradient, while penalized coordinates numerically at
#' zero must satisfy the subgradient condition. The violation is scaled by
#' the number of exceedances and compared with \code{grad_tol}.
#'
#' When \code{lambda = "bic"}, the regularization path is evaluated on a
#' coarse grid followed by a refined grid around the provisional optimum (or
#' on a user-supplied \code{lambda_grid}). The returned model is the path
#' solution at the selected penalty, so the reported fit is exactly the one
#' scored by the BIC.
#'
#' @param df A data frame containing the response and all covariates.
#'   Must have a column named \code{y} (the observations).
#' @param threshold Numeric scalar. The extreme value threshold u.
#' @param loc_vars Character vector of covariate names for the location
#'   parameter. \code{NULL} for stationary location.
#' @param scale_vars Character vector of covariate names for the log-scale
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
#'   or \code{"bic"} to select automatically by BIC along a regularization
#'   path. Ignored when \code{penalty = "none"}.
#' @param penalize_shape Logical. Penalize shape parameter covariates?
#'   Default \code{TRUE}.
#' @param obs_per_year Numeric. Observations per year. E.g. \code{365.25}
#'   for daily data, \code{365.25 / 7} for weekly. Default \code{365.25}.
#' @param active_tol Positive numeric tolerance used for BIC complexity
#'   counting, active-covariate extraction, and default coefficient reporting.
#'   Default \code{1e-2}.
#' @param lambda_scaling Character. Scaling used when \code{lambda = "bic"}:
#'   \code{"common"} (default) applies the same scalar lambda to all parameter
#'   blocks; \code{"gradient"} calibrates the scale and shape penalties
#'   relative to the location-block score magnitude at the initial model.
#'   Ignored for numeric \code{lambda}.
#' @param penalty_weights Coordinate weights multiplying \code{lambda}:
#'   \code{"none"} (default, unit weights), \code{"adaptive"} for the adaptive
#'   LASSO with weights \eqn{1/|\tilde\theta_k|^\gamma} computed from an
#'   initial unpenalized fit (ridge fallback), or a non-negative numeric
#'   vector with one entry per parameter (names, if present, must match the
#'   parameter names). Weights of unpenalized coordinates are ignored.
#' @param adaptive_gamma Positive exponent \eqn{\gamma} of the adaptive
#'   weights. Default \code{1}.
#' @param lambda_grid Optional decreasing numeric vector of base penalties
#'   used for BIC selection instead of the default two-phase grid.
#' @param init Optional numeric starting vector, one entry per parameter.
#' @param grad_tol Positive tolerance for the KKT violation divided by the
#'   number of exceedances. Default \code{5e-4}.
#' @param maxit Integer. Maximum optimizer iterations. Default \code{10000L}.
#' @param calc_hessian Logical. Compute the Hessian of the penalized
#'   objective at the solution? Default \code{FALSE}.
#' @param verbose Logical. Print progress during BIC grid search and
#'   convergence warnings. Default \code{TRUE}.
#'
#' @return An object of class \code{nhpp_fit}. Besides the smooth optimizer
#'   solution (\code{par}) and the operational activity-thresholded solution
#'   (\code{par_oper}), it stores the penalty weights, the scaled KKT
#'   violation (\code{kkt}), and, for \code{lambda = "bic"}, the evaluated
#'   regularization path (\code{path}).
#'
#' @seealso \code{\link{marginalize}}, \code{\link{bootstrap_coef}},
#'   \code{\link{plot.nhpp_fit}}
#'
#' @examples
#' set.seed(42)
#' n <- 2000
#' df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' df$y <- 10 + 2 * df$x1 + rexp(n, 1 / 3)
#'
#' # Unpenalized stationary-shape model
#' fit0 <- fit_nhpp(df, threshold = 18, loc_vars = c("x1", "x2"),
#'                  penalty = "none", verbose = FALSE)
#' coef(fit0)
#'
#' \donttest{
#' # Adaptive LASSO with BIC-selected penalty
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = c("x1", "x2"),
#'                 penalty = "lasso", lambda = "bic",
#'                 penalty_weights = "adaptive", verbose = FALSE)
#' summary(fit)
#' }
#'
#' @export
fit_nhpp <- function(df, threshold,
                     loc_vars        = NULL,
                     scale_vars      = NULL,
                     shape_vars      = NULL,
                     free_vars       = NULL,
                     penalty         = c("none", "lasso", "ridge", "elnet"),
                     alpha           = 0.5,
                     lambda          = "bic",
                     penalize_shape  = TRUE,
                     obs_per_year    = 365.25,
                     active_tol      = 1e-2,
                     lambda_scaling  = c("common", "gradient"),
                     penalty_weights = "none",
                     adaptive_gamma  = 1,
                     lambda_grid     = NULL,
                     init            = NULL,
                     grad_tol        = 5e-4,
                     maxit           = 10000L,
                     calc_hessian    = FALSE,
                     verbose         = TRUE) {

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
  if (!is.numeric(grad_tol) || length(grad_tol) != 1L ||
      !is.finite(grad_tol) || grad_tol <= 0)
    stop("fit_nhpp: `grad_tol` must be a single positive numeric value.")
  if (!is.numeric(adaptive_gamma) || length(adaptive_gamma) != 1L ||
      !is.finite(adaptive_gamma) || adaptive_gamma <= 0)
    stop("fit_nhpp: `adaptive_gamma` must be a single positive numeric value.")
  if (!is.logical(penalize_shape) || length(penalize_shape) != 1L ||
      is.na(penalize_shape))
    stop("fit_nhpp: `penalize_shape` must be TRUE or FALSE.")
  if (!is.logical(calc_hessian) || length(calc_hessian) != 1L ||
      is.na(calc_hessian))
    stop("fit_nhpp: `calc_hessian` must be TRUE or FALSE.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("fit_nhpp: `verbose` must be TRUE or FALSE.")
  if (!is.null(lambda_grid) &&
      (!is.numeric(lambda_grid) || length(lambda_grid) < 2L ||
       any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)))
    stop("fit_nhpp: `lambda_grid` must contain at least two positive finite values.")

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
  par_names <- c(
    paste0("mu.",    colnames(dm$X_mu)),
    paste0("sigma.", colnames(dm$X_sigma)),
    paste0("xi.",    colnames(dm$X_xi))
  )

  if (is.null(init)) {
    y_exc <- y[y > threshold]
    init  <- rep(0, p_mu + p_sig + p_xi)
    init[1L]        <- mean(y_exc)
    init[p_mu + 1L] <- log(max(stats::sd(y_exc), 0.1))
  } else {
    if (!is.numeric(init) || length(init) != length(par_names) ||
        any(!is.finite(init)))
      stop(sprintf("fit_nhpp: `init` must be a finite numeric vector of length %d.",
                   length(par_names)))
    if (!is.null(names(init))) {
      if (!setequal(names(init), par_names))
        stop("fit_nhpp: names of `init` must match the model parameter names.")
      init <- init[par_names]
    }
  }
  init <- as.numeric(init)
  names(init) <- par_names

  # ---- penalty weights -------------------------------------------------------
  weights_source <- "none"
  pen_weights    <- NULL
  if (is.character(penalty_weights)) {
    if (length(penalty_weights) != 1L ||
        !penalty_weights %in% c("none", "adaptive"))
      stop("fit_nhpp: `penalty_weights` must be \"none\", \"adaptive\", or a numeric vector.")
    if (penalty_weights == "adaptive" && penalty != "none") {
      aw <- .adaptive_weights(dm, y, threshold, init, penalize_shape,
                              obs_per_year, gamma = adaptive_gamma,
                              maxit = maxit, grad_tol = grad_tol,
                              active_tol = active_tol)
      pen_weights    <- aw$weights
      weights_source <- paste0("adaptive: ", aw$source)
      if (verbose)
        message("fit_nhpp: adaptive weights from ", aw$source, " fit.")
    }
  } else {
    if (!is.numeric(penalty_weights) ||
        length(penalty_weights) != length(par_names) ||
        any(!is.finite(penalty_weights)) || any(penalty_weights < 0))
      stop(sprintf(
        "fit_nhpp: numeric `penalty_weights` must be %d finite non-negative values.",
        length(par_names)))
    if (!is.null(names(penalty_weights))) {
      if (!setequal(names(penalty_weights), par_names))
        stop("fit_nhpp: names of `penalty_weights` must match the parameter names.")
      penalty_weights <- penalty_weights[par_names]
    }
    pen_weights <- stats::setNames(as.numeric(penalty_weights), par_names)
    weights_source <- "user"
  }

  # ---- lambda selection --------------------------------------------------------
  path <- NULL
  start <- init
  if (is.null(lambda_resolved)) {
    if (verbose)
      message("fit_nhpp: running BIC lambda selection (penalty = '",
              penalty, "')...")
    sel <- .select_lambda_bic(
      dm              = dm,
      y               = y,
      threshold       = threshold,
      alpha           = alpha,
      penalize_shape  = penalize_shape,
      init            = init,
      obs_per_year    = obs_per_year,
      active_tol      = active_tol,
      lambda_scaling  = lambda_scaling,
      maxit           = maxit,
      verbose         = verbose,
      penalty_weights = pen_weights,
      lambda_grid     = lambda_grid,
      grad_tol        = grad_tol
    )
    lambda_resolved <- sel$lambda
    start <- sel$par
    path  <- sel$path
    if (verbose)
      message(sprintf("fit_nhpp: selected lambda = %s",
                      paste(sprintf("%s=%.5g", names(lambda_resolved),
                                    lambda_resolved), collapse = ", ")))
  } else if (any(lambda_resolved > 0)) {
    start <- .warmstart(init, dm, y, threshold, lambda_resolved, alpha,
                        penalize_shape, obs_per_year = obs_per_year,
                        penalty_weights = pen_weights)
  }

  res <- .fit_at_lambda(dm, y, threshold,
                        lambda          = lambda_resolved,
                        alpha           = alpha,
                        penalize_shape  = penalize_shape,
                        init            = start,
                        maxit           = maxit,
                        calc_hessian    = calc_hessian,
                        obs_per_year    = obs_per_year,
                        penalty_weights = pen_weights,
                        grad_tol        = grad_tol,
                        active_tol      = active_tol)

  if (!res$converged && verbose)
    warning("fit_nhpp: optimizer did not converge (KKT check failed).")

  par_hat <- res$par
  names(par_hat) <- par_names

  pen_idx_all <- .penalized_index(dm, penalize_shape)
  par_oper <- par_hat
  small <- pen_idx_all[abs(par_oper[pen_idx_all]) < active_tol]
  if (length(small) > 0L) par_oper[small] <- 0
  nllh_oper <- pp_nllh(par_oper, dm, y, threshold,
                       lambda = 0, alpha = alpha,
                       pen_xi = penalize_shape,
                       obs_per_year = obs_per_year)

  fitted_from <- function(p) list(
    mu    = as.numeric(dm$X_mu    %*% p[seq_len(p_mu)]),
    sigma = exp(as.numeric(dm$X_sigma %*% p[p_mu + seq_len(p_sig)])),
    xi    = as.numeric(dm$X_xi    %*% p[p_mu + p_sig + seq_len(p_xi)])
  )

  new_nhpp_fit(
    par             = par_hat,
    par_oper        = par_oper,
    dm              = dm,
    threshold       = threshold,
    nllh_pen        = res$nllh_pen,
    nllh_raw        = res$nllh_raw,
    nllh_oper       = nllh_oper,
    lambda          = lambda_resolved,
    alpha           = alpha,
    penalty         = penalty,
    penalize_shape  = penalize_shape,
    lambda_scaling  = lambda_scaling,
    hessian         = res$hessian,
    fitted          = fitted_from(par_hat),
    fitted_oper     = fitted_from(par_oper),
    obs_per_year    = obs_per_year,
    active_tol      = active_tol,
    converged       = res$converged,
    n_exc           = n_exc,
    penalty_weights = pen_weights,
    weights_source  = weights_source,
    kkt             = res$kkt,
    path            = path,
    y               = y,
    grad_tol        = grad_tol,
    adaptive_gamma  = adaptive_gamma
  )
}

# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

.format_lambda <- function(lambda) {
  if (length(lambda) == 1L) return(sprintf("%.5g", lambda))
  nm <- if (is.null(names(lambda))) c("mu", "sigma", "xi") else names(lambda)
  paste(sprintf("%s=%.5g", nm, lambda), collapse = ", ")
}

#' Print a fitted nhpp_fit model
#'
#' @param x An \code{nhpp_fit} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' print(fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                verbose = FALSE))
#'
#' @export
print.nhpp_fit <- function(x, ...) {
  cat("-- nhpp_fit --------------------------------------\n")
  cat(sprintf("  Threshold    : %.4g\n",  x$threshold))
  cat(sprintf("  Penalty      : %s\n",    x$penalty))
  cat(sprintf("  Lambda       : %s\n",    .format_lambda(x$lambda)))
  if (!is.null(x$weights_source) && !identical(x$weights_source, "none"))
    cat(sprintf("  Weights      : %s\n",  x$weights_source))
  cat(sprintf("  Alpha        : %s\n",    .format_lambda(x$alpha)))
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

#' Extract coefficients of a fitted nhpp_fit model
#'
#' @param object An \code{nhpp_fit} object.
#' @param operational Logical. If \code{TRUE}, return the operational
#'   estimator, in which penalized coefficients below the activity tolerance
#'   are set to zero. Default \code{FALSE} (smooth optimizer solution).
#' @param ... Ignored.
#'
#' @return Named numeric vector of coefficients.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' coef(fit)
#'
#' @export
coef.nhpp_fit <- function(object, operational = FALSE, ...) {
  if (isTRUE(operational) && !is.null(object$par_oper)) object$par_oper else
    object$par
}
