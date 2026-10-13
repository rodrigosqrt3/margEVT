# =============================================================================
# bootstrap.R
# Parametric bootstrap confidence intervals for return levels and coefficients.
#
# Two user-facing functions:
#   - bootstrap_rl()   : CI for return levels
#   - bootstrap_coef() : CI for model coefficients
#
# Both refit at the SAME lambda, penalty weights, and activity tolerance as
# the original fit, starting from the original estimate. The resulting
# intervals are conditional on the selected specification; they do not
# propagate model-selection uncertainty.
#
# Every replicate b draws its random numbers from its own L'Ecuyer-CMRG
# stream, so the results are identical for any number of cores.
#
# Internal workhorse:
#   - .simulate_exceedances() : simulate new exceedances from fitted GPD
#   - .refit_boot()           : refit model on one bootstrap sample
#   - .boot_replicates()      : run the replicates (serially or in parallel)
# =============================================================================

# -----------------------------------------------------------------------------
# .simulate_exceedances()
# Given fitted parameters and exceedance indices, simulate new exceedances
# from the implied GPD via quantile inversion.
# Internal — not exported.
# -----------------------------------------------------------------------------
# GPD scale at the threshold, sigma + xi * (u - mu), at the exceedance times.
.gpd_scale_at_threshold <- function(mu_t, sigma_t, xi_t, threshold, exc_idx) {
  sig_gpd <- sigma_t[exc_idx] + xi_t[exc_idx] * (threshold - mu_t[exc_idx])
  if (any(!is.finite(sig_gpd)) || any(sig_gpd <= 0))
    stop("bootstrap: invalid GPD scale at threshold - check model parameters.")
  sig_gpd
}

.simulate_exceedances <- function(mu_t, sigma_t, xi_t, threshold, exc_idx) {

  xi_e    <- xi_t[exc_idx]
  sig_gpd <- .gpd_scale_at_threshold(mu_t, sigma_t, xi_t, threshold, exc_idx)

  u_unif <- pmin(pmax(stats::runif(length(exc_idx)),
                      .Machine$double.eps),
                 1 - .Machine$double.eps)
  out <- numeric(length(exc_idx))
  gumbel <- abs(xi_e) < 1e-6
  out[gumbel] <- -sig_gpd[gumbel] * log1p(-u_unif[gumbel])
  out[!gumbel] <- sig_gpd[!gumbel] *
    expm1(-xi_e[!gumbel] * log1p(-u_unif[!gumbel])) / xi_e[!gumbel]
  if (any(!is.finite(out)) || any(out < 0))
    stop("bootstrap: generated an invalid GPD excess.")
  out
}

# -----------------------------------------------------------------------------
# .refit_boot()
# Refit the model on new data with the original design and penalty settings.
#   init            : starting values (default: the original estimate)
#   lambda          : penalty level (default: the original one; "bic" re-selects)
#   penalty_weights : NULL keeps the original weights; "adaptive" recomputes
#                     them on `df_boot`.
# Returns an nhpp_fit or NULL on failure.
# Internal — not exported.
# -----------------------------------------------------------------------------
.refit_boot <- function(df_boot, fit, init = fit$par, lambda = fit$lambda,
                        penalty_weights = NULL) {

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

  if (is.null(penalty_weights))
    penalty_weights <- if (is.null(fit$penalty_weights)) "none" else
      fit$penalty_weights

  tryCatch(
    fit_nhpp(df_boot, fit$threshold,
             loc_vars        = loc_vars,
             scale_vars      = sig_vars,
             shape_vars      = xi_vars,
             free_vars       = free_vars,
             penalty         = fit$penalty,
             alpha           = fit$alpha,
             lambda          = lambda,
             penalize_shape  = fit$penalize_shape,
             obs_per_year    = fit$obs_per_year,
             active_tol      = if (is.null(fit$active_tol)) 1e-2 else fit$active_tol,
             lambda_scaling  = if (is.null(fit$lambda_scaling)) "common" else
               fit$lambda_scaling,
             penalty_weights = penalty_weights,
             adaptive_gamma  = if (is.null(fit$adaptive_gamma)) 1 else
               fit$adaptive_gamma,
             init            = init,
             grad_tol        = if (is.null(fit$grad_tol)) 5e-4 else fit$grad_tol,
             maxit           = 5000L,
             calc_hessian    = FALSE,
             verbose         = FALSE),
    error = function(e) NULL
  )
}

# -----------------------------------------------------------------------------
# .boot_replicates()
# Runs R parametric bootstrap replicates. `stat(fit_b, df_boot)` returns a
# numeric vector of length `n_stat` or NULL; the result is an R x n_stat
# matrix with NA rows for failed replicates.
# Internal — not exported.
# -----------------------------------------------------------------------------
.boot_replicates <- function(fit, data, R, seed, cores, stat, n_stat,
                             fun_name) {
  thr     <- fit$threshold
  exc_idx <- which(data$y > thr)
  if (length(exc_idx) < 5L)
    stop(sprintf("%s: fewer than 5 exceedances - bootstrap unreliable.",
                 fun_name))

  params  <- predict_params(fit)
  # Fail early, in the calling process, when the GPD scale is invalid.
  .gpd_scale_at_threshold(params$mu, params$sigma, params$xi, thr, exc_idx)

  streams <- .rng_streams(R, seed)
  one_replicate <- function(stream) .with_rng_stream(stream, {
    excess_sim         <- .simulate_exceedances(params$mu, params$sigma,
                                                params$xi, thr, exc_idx)
    df_boot            <- data
    df_boot$y          <- thr - stats::runif(nrow(data), 0.1, 5)
    df_boot$y[exc_idx] <- thr + excess_sim
    fit_b <- .refit_boot(df_boot, fit)
    if (is.null(fit_b) || !isTRUE(fit_b$converged)) return(NULL)
    stat(fit_b, df_boot)
  })

  out <- .par_lapply(streams, one_replicate, cores = cores)
  mat <- matrix(NA_real_, nrow = R, ncol = n_stat)
  for (b in seq_len(R))
    if (!is.null(out[[b]])) mat[b, ] <- out[[b]]
  mat
}

.validate_boot_inputs <- function(fit, data, R, level, seed, fun) {
  if (!inherits(fit, "nhpp_fit"))
    stop(sprintf("%s: `fit` must be an nhpp_fit object.", fun))
  if (!is.data.frame(data))
    stop(sprintf("%s: `data` must be a data frame.", fun))
  if (!"y" %in% names(data) || !is.numeric(data$y) ||
      any(!is.finite(data$y)))
    stop(sprintf("%s: `data$y` must contain only finite numeric values.", fun))
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 1)
    stop(sprintf("%s: `R` must be a positive integer.", fun))
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop(sprintf("%s: `level` must lie strictly between 0 and 1.", fun))
  .validate_seed(seed, fun)
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# bootstrap_rl()
# -----------------------------------------------------------------------------

#' Parametric bootstrap confidence intervals for return levels
#'
#' Simulates \code{R} bootstrap samples from the fitted model, refits at the
#' original model design, penalty, lambda, penalty weights, and activity
#' tolerance (starting from the original estimate), recomputes return levels
#' via \code{\link{marginalize}}, and returns empirical percentile intervals.
#' Exceedance count and times are held fixed, so the result is a nominal
#' conditional summary rather than unconditional or post-selection
#' inference.
#'
#' Each replicate uses its own L'Ecuyer-CMRG random-number stream derived
#' from \code{seed}, so the output does not depend on \code{cores}. The
#' caller's random-number state is left untouched when \code{seed} is
#' supplied.
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
#'   \code{scenarios}, \code{fixed}, \code{interactions}).
#' @param level Numeric. Confidence level. Default \code{0.95}.
#' @param seed Seed of the replicate streams, or \code{NULL} to derive it
#'   from the current random-number state. Default \code{42L}.
#' @param cores Positive integer number of worker processes. Default \code{1L}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A data frame with columns \code{TR}, \code{RL_est} (the return
#'   level from the original fitted model under the requested approach),
#'   \code{CI_low}, \code{CI_high}, \code{n_ok} (number of successful
#'   replicates). The \code{R x length(TRs)} matrix of replicate return
#'   levels is attached as attribute \code{"replicates"}.
#'
#' @seealso \code{\link{bootstrap_coef}}, \code{\link{marginalize}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' \donttest{
#' bootstrap_rl(fit, df, TRs = c(10, 50), R = 50, approach = "A",
#'              verbose = FALSE)
#' }
#'
#' @export
bootstrap_rl <- function(fit, data,
                         TRs              = c(2, 5, 10, 20, 50, 100),
                         R                = 200L,
                         approach         = "C",
                         marginalize_args = list(),
                         level            = 0.95,
                         seed             = 42L,
                         cores            = 1L,
                         verbose          = TRUE) {

  .validate_boot_inputs(fit, data, R, level, seed, "bootstrap_rl")
  if (!is.character(approach) || length(approach) != 1L ||
      !approach %in% c("A", "B", "C"))
    stop("bootstrap_rl: `approach` must be one of 'A', 'B', 'C'.")
  if (!is.numeric(TRs) || length(TRs) < 1L ||
      any(!is.finite(TRs)) || any(TRs <= 1))
    stop("bootstrap_rl: `TRs` must contain finite return periods greater than 1.")
  cores <- .validate_cores(cores, "bootstrap_rl")
  R <- as.integer(R)

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

  rl_of <- function(f, d) {
    args <- utils::modifyList(
      effective_marginalize_args,
      list(fit = f, data = d, TRs = TRs, approaches = approach)
    )
    res <- tryCatch(do.call(marginalize, args), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    vapply(TRs, function(tr) {
      value <- res$RL[res$TR == tr]
      if (length(value) == 0L) NA_real_ else value[1L]
    }, numeric(1L))
  }

  point_est <- rl_of(fit, data)
  if (is.null(point_est)) point_est <- rep(NA_real_, length(TRs))

  if (verbose)
    message(sprintf("bootstrap_rl: R=%d | approach=%s | level=%.0f%% | cores=%d",
                    R, approach, level * 100, cores))

  boot_mat <- .boot_replicates(fit, data, R, seed, cores, rl_of,
                               length(TRs), "bootstrap_rl")
  colnames(boot_mat) <- paste0("T", TRs)

  if (verbose)
    message(sprintf("  %d / %d replicates succeeded",
                    sum(stats::complete.cases(boot_mat)), R))

  alpha_tail <- (1 - level) / 2
  out <- do.call(rbind, lapply(seq_along(TRs), function(ti) {
    v <- stats::na.omit(boot_mat[, ti])
    data.frame(
      TR       = TRs[ti],
      RL_est   = round(point_est[ti], 4L),
      CI_low   = if (length(v) > 0L)
        round(unname(stats::quantile(v, alpha_tail)), 4L) else NA_real_,
      CI_high  = if (length(v) > 0L)
        round(unname(stats::quantile(v, 1 - alpha_tail)), 4L) else NA_real_,
      n_ok     = length(v)
    )
  }))
  attr(out, "replicates") <- boot_mat
  out
}


# -----------------------------------------------------------------------------
# bootstrap_coef()
# -----------------------------------------------------------------------------

#' Parametric bootstrap confidence intervals for model coefficients
#'
#' Simulates \code{R} bootstrap samples, refits at the original penalty,
#' lambda, and penalty weights, and returns empirical percentile intervals
#' for each operational coefficient. See \code{\link{bootstrap_rl}} for the
#' resampling scheme and the random-number streams.
#'
#' @inheritParams bootstrap_rl
#'
#' @return A data frame with columns \code{parameter}, \code{estimate},
#'   \code{CI_low}, \code{CI_high}, \code{n_ok}, with one row for every
#'   coefficient in the fitted model. Intervals require at least five
#'   successful replicates. The replicate matrix is attached as attribute
#'   \code{"replicates"}.
#'
#' @seealso \code{\link{bootstrap_rl}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' \donttest{
#' bootstrap_coef(fit, df, R = 50, verbose = FALSE)
#' }
#'
#' @export
bootstrap_coef <- function(fit, data,
                           R       = 200L,
                           level   = 0.95,
                           seed    = 42L,
                           cores   = 1L,
                           verbose = TRUE) {

  .validate_boot_inputs(fit, data, R, level, seed, "bootstrap_coef")
  cores <- .validate_cores(cores, "bootstrap_coef")
  R <- as.integer(R)
  par_names <- names(fit$par)

  if (verbose)
    message(sprintf("bootstrap_coef: R=%d | level=%.0f%% | cores=%d",
                    R, level * 100, cores))

  coef_of <- function(f, d) {
    par_b <- if (is.null(f$par_oper)) f$par else f$par_oper
    par_b[par_names]
  }
  boot_coef <- .boot_replicates(fit, data, R, seed, cores, coef_of,
                                length(par_names), "bootstrap_coef")
  colnames(boot_coef) <- par_names

  if (verbose)
    message(sprintf("  %d / %d replicates succeeded",
                    sum(stats::complete.cases(boot_coef)), R))

  alpha_tail <- (1 - level) / 2
  est <- if (is.null(fit$par_oper)) fit$par else fit$par_oper
  out <- do.call(rbind, lapply(par_names, function(nm) {
    v <- stats::na.omit(boot_coef[, nm])
    data.frame(
      parameter = nm,
      estimate  = round(unname(est[nm]), 5L),
      CI_low    = if (length(v) >= 5L)
        round(unname(stats::quantile(v, alpha_tail)), 5L) else NA_real_,
      CI_high   = if (length(v) >= 5L)
        round(unname(stats::quantile(v, 1 - alpha_tail)), 5L) else NA_real_,
      n_ok      = length(v),
      row.names = NULL
    )
  }))
  attr(out, "replicates") <- boot_coef
  out
}
