# =============================================================================
# bootstrap.R
# Parametric bootstrap confidence intervals for return levels and coefficients.
#
# Two user-facing functions:
#   - bootstrap_rl()   : CI for return levels
#   - bootstrap_coef() : CI for model coefficients
#
# Both refit at the SAME lambda/alpha and activity tolerance as the original
# fit. The resulting intervals are conditional on the selected specification;
# they do not propagate model-selection uncertainty.
#
# Internal workhorse:
#   - .simulate_exceedances() : simulate new exceedances from fitted GPD
#   - .refit_boot()           : refit model on one bootstrap sample
# =============================================================================

# -----------------------------------------------------------------------------
# .simulate_exceedances()
# Given fitted parameters and exceedance indices, simulate new exceedances
# from the implied GPD via quantile inversion.
# Internal — not exported.
# -----------------------------------------------------------------------------
.simulate_exceedances <- function(mu_t, sigma_t, xi_t, threshold, exc_idx) {

  mu_e    <- mu_t[exc_idx]
  sig_e   <- sigma_t[exc_idx]
  xi_e    <- xi_t[exc_idx]

  # GPD scale at threshold: sigma_gpd = sigma + xi * (u - mu)
  sig_gpd <- sig_e + xi_e * (threshold - mu_e)

  if (any(!is.finite(sig_gpd)) || any(sig_gpd <= 0))
    stop("bootstrap: invalid GPD scale at threshold - check model parameters.")

  u_unif <- stats::runif(length(exc_idx))

  ifelse(
    abs(xi_e) < 1e-6,
    -sig_gpd * log(u_unif),                              # Exponential limit
    (sig_gpd / xi_e) * (u_unif^(-xi_e) - 1)             # GPD quantile
  )
}

# -----------------------------------------------------------------------------
# .refit_boot()
# Refit the model on one bootstrap sample at the original lambda/alpha.
# Returns an nhpp_fit or NULL on failure.
# Internal — not exported.
# -----------------------------------------------------------------------------
.refit_boot <- function(df_boot, fit) {

  dm       <- fit$dm
  loc_vars  <- colnames(dm$X_mu)[colnames(dm$X_mu) != "(Intercept)"]
  sig_vars  <- colnames(dm$X_sigma)[colnames(dm$X_sigma) != "(Intercept)"]
  xi_vars   <- colnames(dm$X_xi)[colnames(dm$X_xi)   != "(Intercept)"]

  loc_vars  <- if (length(loc_vars) == 0L) NULL else loc_vars
  sig_vars  <- if (length(sig_vars) == 0L) NULL else sig_vars
  xi_vars   <- if (length(xi_vars)  == 0L) NULL else xi_vars

  # Reconstruct the union of unpenalized columns across all parameter blocks.
  # This preserves mandatory seasonal terms and any user-specified free
  # variables in exactly the same model design used by the original fit.
  free_in_block <- function(X, idx_pen) {
    cols <- colnames(X)
    setdiff(cols, c("(Intercept)", cols[idx_pen]))
  }
  free_vars <- unique(c(
    free_in_block(dm$X_mu, dm$idx_pen_mu),
    free_in_block(dm$X_sigma, dm$idx_pen_sigma),
    free_in_block(dm$X_xi, dm$idx_pen_xi)
  ))
  free_vars    <- if (length(free_vars) == 0L) NULL else free_vars

  tryCatch(
    fit_nhpp(df_boot, fit$threshold,
             loc_vars       = loc_vars,
             scale_vars     = sig_vars,
             shape_vars     = xi_vars,
             free_vars      = free_vars,
             penalty        = fit$penalty,
             alpha          = fit$alpha,
             lambda         = fit$lambda,
             penalize_shape = fit$penalize_shape,
             obs_per_year   = fit$obs_per_year,
             active_tol     = if (is.null(fit$active_tol)) 1e-2 else fit$active_tol,
             maxit          = 5000L,
             calc_hessian   = FALSE,
             verbose        = FALSE),
    error = function(e) NULL
  )
}


# -----------------------------------------------------------------------------
# bootstrap_rl()
# -----------------------------------------------------------------------------

#' Parametric bootstrap confidence intervals for return levels
#'
#' Simulates \code{R} bootstrap samples from the fitted model, refits at the
#' original model design, penalty, lambda, and activity tolerance, recomputes
#' return levels via \code{\link{marginalize}}, and returns empirical quantile
#' intervals. Exceedance count and times are held fixed, so the result is a
#' nominal conditional summary rather than unconditional or post-selection
#' inference.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data The data frame used to fit \code{fit}. Must contain column
#'   \code{y} and all covariates.
#' @param TRs Numeric vector of return periods. Default \code{c(2,5,10,20,50,100)}.
#' @param R Integer. Number of bootstrap replicates. Default \code{200L}.
#' @param approach Character. Which marginalization approach to use for each
#'   replicate: \code{"A"}, \code{"B"}, or \code{"C"}. Default \code{"C"}.
#' @param marginalize_args Named list of additional arguments passed to
#'   \code{\link{marginalize}} (e.g. \code{n_boot}, \code{seed},
#'   \code{scenarios}, \code{interactions}).
#' @param level Numeric. Confidence level. Default \code{0.95}.
#' @param seed Integer. Random seed. Default \code{42L}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A data frame with columns \code{TR}, \code{RL_est} (the return
#'   level from the original fitted model under the requested approach),
#'   \code{CI_low}, \code{CI_high}, \code{n_ok} (number of successful
#'   replicates).
#'
#' @export
bootstrap_rl <- function(fit, data,
                         TRs              = c(2, 5, 10, 20, 50, 100),
                         R                = 200L,
                         approach         = "C",
                         marginalize_args = list(),
                         level            = 0.95,
                         seed             = 42L,
                         verbose          = TRUE) {

  if (!inherits(fit, "nhpp_fit"))
    stop("bootstrap_rl: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("bootstrap_rl: `data` must be a data frame.")
  if (!approach %in% c("A", "B", "C"))
    stop("bootstrap_rl: `approach` must be one of 'A', 'B', 'C'.")
  if (!is.numeric(TRs) || length(TRs) < 1L ||
      any(!is.finite(TRs)) || any(TRs <= 1))
    stop("bootstrap_rl: `TRs` must contain finite return periods greater than 1.")
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 1)
    stop("bootstrap_rl: `R` must be a positive integer.")
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("bootstrap_rl: `level` must lie strictly between 0 and 1.")
  R <- as.integer(R)

  set.seed(seed)

  y       <- data$y
  thr     <- fit$threshold
  exc_idx <- which(y > thr)

  if (length(exc_idx) < 5L)
    stop("bootstrap_rl: fewer than 5 exceedances - bootstrap unreliable.")

  params  <- predict_params(fit)
  mu_t    <- params$mu
  sigma_t <- params$sigma
  xi_t    <- params$xi

  alpha_tail <- (1 - level) / 2
  boot_mat   <- matrix(NA_real_, nrow = R, ncol = length(TRs),
                       dimnames = list(NULL, paste0("T", TRs)))

  effective_marginalize_args <- marginalize_args
  if (approach == "A" && !"scenarios" %in% names(marginalize_args)) {
    ac <- active_covariates(fit)
    interaction_inputs <- if ("interactions" %in% names(marginalize_args))
      unique(unlist(marginalize_args$interactions, use.names = FALSE)) else
        character(0L)
    scenario_vars <- unique(c(ac, interaction_inputs))
    scenario_vars <- scenario_vars[scenario_vars %in% names(data)]
    default_sc <- stats::setNames(
      lapply(scenario_vars, function(v) mean(data[[v]], na.rm = TRUE)),
      scenario_vars
    )
    effective_marginalize_args$scenarios <- list(mean = default_sc)
  }

  point_args <- utils::modifyList(
    effective_marginalize_args,
    list(fit = fit, data = data, TRs = TRs, approaches = approach)
  )
  point_result <- tryCatch(do.call(marginalize, point_args),
                           error = function(e) NULL)
  point_est <- vapply(TRs, function(tr) {
    if (is.null(point_result)) return(NA_real_)
    value <- point_result$RL[point_result$TR == tr]
    if (length(value) == 0L) NA_real_ else value[1L]
  }, numeric(1L))

  if (verbose)
    message(sprintf("bootstrap_rl: R=%d | approach=%s | level=%.0f%%",
                    R, approach, level * 100))

  for (b in seq_len(R)) {

    # Simulate new exceedances from the fitted GPD
    excess_sim <- .simulate_exceedances(mu_t, sigma_t, xi_t, thr, exc_idx)

    df_boot              <- data
    df_boot$y            <- thr - stats::runif(nrow(data), 0.1, 5)
    df_boot$y[exc_idx]   <- thr + excess_sim

    fit_b <- .refit_boot(df_boot, fit)
    if (is.null(fit_b) || !fit_b$converged) next

    marg_args <- utils::modifyList(
      effective_marginalize_args,
      list(fit = fit_b, data = df_boot, TRs = TRs, approaches = approach)
    )
    rl_b <- tryCatch(do.call(marginalize, marg_args), error = function(e) NULL)
    if (is.null(rl_b)) next

    for (ti in seq_along(TRs))
      boot_mat[b, ti] <- rl_b$RL[rl_b$TR == TRs[ti]][1L]

    if (verbose && (b %% 25L == 0L || b == R))
      message(sprintf("  replicate %d / %d", b, R))
  }

  do.call(rbind, lapply(seq_along(TRs), function(ti) {
    v <- stats::na.omit(boot_mat[, ti])
    data.frame(
      TR       = TRs[ti],
      RL_est   = round(point_est[ti], 4L),
      CI_low   = if (length(v) > 0L) round(stats::quantile(v, alpha_tail),  4L) else NA_real_,
      CI_high  = if (length(v) > 0L) round(stats::quantile(v, 1 - alpha_tail), 4L) else NA_real_,
      n_ok     = length(v)
    )
  }))
}


# -----------------------------------------------------------------------------
# bootstrap_coef()
# -----------------------------------------------------------------------------

#' Parametric bootstrap confidence intervals for model coefficients
#'
#' Simulates \code{R} bootstrap samples, refits at the original penalty and
#' lambda, and returns empirical quantile CIs for each coefficient.
#'
#' @inheritParams bootstrap_rl
#'
#' @return A data frame with columns \code{parameter}, \code{estimate},
#'   \code{CI_low}, \code{CI_high}, \code{n_ok}, with one row for every
#'   coefficient in the fitted model.
#'
#' @export
bootstrap_coef <- function(fit, data,
                           R       = 200L,
                           level   = 0.95,
                           seed    = 42L,
                           verbose = TRUE) {

  if (!inherits(fit, "nhpp_fit"))
    stop("bootstrap_coef: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("bootstrap_coef: `data` must be a data frame.")
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 1)
    stop("bootstrap_coef: `R` must be a positive integer.")
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("bootstrap_coef: `level` must lie strictly between 0 and 1.")
  R <- as.integer(R)

  set.seed(seed)

  y       <- data$y
  thr     <- fit$threshold
  exc_idx <- which(y > thr)

  if (length(exc_idx) < 5L)
    stop("bootstrap_coef: fewer than 5 exceedances - bootstrap unreliable.")

  params  <- predict_params(fit)
  mu_t    <- params$mu
  sigma_t <- params$sigma
  xi_t    <- params$xi

  alpha_tail <- (1 - level) / 2
  boot_coef  <- matrix(NA_real_, nrow = R, ncol = length(fit$par),
                       dimnames = list(NULL, names(fit$par)))

  if (verbose)
    message(sprintf("bootstrap_coef: R=%d | level=%.0f%%", R, level * 100))

  for (b in seq_len(R)) {

    excess_sim <- .simulate_exceedances(mu_t, sigma_t, xi_t, thr, exc_idx)

    df_boot            <- data
    df_boot$y          <- thr - stats::runif(nrow(data), 0.1, 5)
    df_boot$y[exc_idx] <- thr + excess_sim

    fit_b <- .refit_boot(df_boot, fit)
    if (is.null(fit_b) || !fit_b$converged) next

    common <- intersect(names(fit$par), names(fit_b$par))
    boot_coef[b, common] <- fit_b$par[common]

    if (verbose && (b %% 25L == 0L || b == R))
      message(sprintf("  replicate %d / %d", b, R))
  }

  do.call(rbind, lapply(names(fit$par), function(nm) {
    v <- stats::na.omit(boot_coef[, nm])
    data.frame(
      parameter = nm,
      estimate  = round(fit$par[nm],                        5L),
      CI_low    = if (length(v) >= 5L)
        round(stats::quantile(v, alpha_tail),       5L) else NA_real_,
      CI_high   = if (length(v) >= 5L)
        round(stats::quantile(v, 1 - alpha_tail),   5L) else NA_real_,
      n_ok      = length(v),
      row.names = NULL
    )
  }))
}
