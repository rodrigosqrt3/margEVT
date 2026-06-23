# =============================================================================
# bootstrap.R
# Parametric bootstrap confidence intervals for return levels and coefficients.
#
# Two user-facing functions:
#   - bootstrap_rl()   : CI for return levels
#   - bootstrap_coef() : CI for model coefficients
#
# Both refit at the SAME lambda/alpha as the original fit — never re-run
# BIC selection on bootstrap samples (that would be double shrinkage and
# produce artificially narrow CIs).
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
  loc_vars  <- colnames(dm$X_mu)[!colnames(dm$X_mu) %in%
                                   c("(Intercept)", "cos1", "sen1", "cos2", "sen2")]
  sig_vars  <- colnames(dm$X_sigma)[colnames(dm$X_sigma) != "(Intercept)"]
  xi_vars   <- colnames(dm$X_xi)[colnames(dm$X_xi)   != "(Intercept)"]

  loc_vars  <- if (length(loc_vars) == 0L) NULL else loc_vars
  sig_vars  <- if (length(sig_vars) == 0L) NULL else sig_vars
  xi_vars   <- if (length(xi_vars)  == 0L) NULL else xi_vars

  # free_vars: columns in X_mu that are NOT in idx_pen_mu
  # (i.e. they were declared free in the original fit)
  all_mu_cols  <- colnames(dm$X_mu)
  penalized_mu <- all_mu_cols[dm$idx_pen_mu]
  free_vars    <- setdiff(all_mu_cols,
                          c("(Intercept)", penalized_mu))
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
#' original penalty and lambda, recomputes return levels via
#' \code{\link{marginalize}}, and returns empirical quantile CIs.
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
#' @return A data frame with columns \code{TR}, \code{RL_est},
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

    # For approach A with no user-supplied scenarios: build a default
    # scenario from df_boot column means to avoid spurious warnings
    marg_args_b <- marginalize_args
    if (approach == "A" && !"scenarios" %in% names(marginalize_args)) {
      ac <- active_covariates(fit_b)
      if (!is.null(ac)) {
        ac_in_boot <- ac[ac %in% names(df_boot)]
        names(ac_in_boot) <- ac_in_boot
        default_sc <- lapply(ac_in_boot, function(v)
          mean(df_boot[[v]], na.rm = TRUE))
        marg_args_b <- c(marg_args_b,
                         list(scenarios = list(mean = default_sc)))
      }
    }

    marg_args <- c(
      list(fit       = fit_b,
           data      = df_boot,
           TRs       = TRs,
           approaches = approach),
      marg_args_b
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
      RL_est   = if (length(v) > 0L) round(mean(v),                    4L) else NA_real_,
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
#'   \code{CI_low}, \code{CI_high}, \code{n_ok}. Only parameters with
#'   non-negligible estimates (\code{|estimate| > 1e-4}) are returned.
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
