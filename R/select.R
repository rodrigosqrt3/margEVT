# =============================================================================
# select.R
# BIC-based two-phase lambda selection for penalized NHPP models.
#
# .select_lambda_bic() is internal — called by fit_nhpp() when lambda = "bic".
# Users never call this directly.
#
# Algorithm:
#   Phase 1: coarse log-spaced grid, early stopping when BIC valley is found.
#   Phase 2: fine grid bracketing the coarse optimum.
#   Block-calibrated lambda: separate lambda per parameter block (mu, sigma, xi)
#   scaled by the gradient magnitude at the null model, so each block is
#   penalized proportionally to its natural scale.
# =============================================================================

#' BIC-based lambda selection for penalized NHPP models (internal)
#'
#' @param dm Design matrix list from \code{build_design_matrices()}.
#' @param y Numeric response vector.
#' @param threshold Numeric scalar threshold.
#' @param alpha Elastic-net mixing parameter.
#' @param penalize_shape Logical.
#' @param init Named numeric starting parameter vector.
#' @param obs_per_year Numeric.
#' @param maxit Integer.
#' @param verbose Logical.
#'
#' @return Named numeric vector \code{c(mu=, sigma=, xi=)} of selected lambdas.
#'
#' @keywords internal
.select_lambda_bic <- function(dm, y, threshold, alpha,
                               penalize_shape, init,
                               obs_per_year = 365.25,
                               maxit        = 5000L,
                               verbose      = TRUE) {

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)
  n_exc <- sum(y > threshold, na.rm = TRUE)

  pen_idx_all <- c(
    dm$idx_pen_mu,
    p_mu + dm$idx_pen_sigma,
    if (penalize_shape) p_mu + p_sig + dm$idx_pen_xi else integer(0L)
  )

  # Block-calibrated lambda ratios
  # Scale each block's lambda by its gradient magnitude at the null model,
  # so that all blocks are penalized proportionally to their natural scale.
  null_grad <- pp_grad(init, dm, y, threshold,
                       lambda = 0, alpha = 1,
                       pen_xi = penalize_shape,
                       obs_per_year = obs_per_year)

  g_mu    <- null_grad[dm$idx_pen_mu]
  g_sigma <- null_grad[p_mu + dm$idx_pen_sigma]
  g_xi    <- if (length(dm$idx_pen_xi) > 0L)
    null_grad[p_mu + p_sig + dm$idx_pen_xi] else numeric(0L)

  scale_mu    <- sqrt(mean(g_mu^2))
  scale_sigma <- if (length(g_sigma) > 0L) sqrt(mean(g_sigma^2)) else scale_mu
  scale_xi    <- if (length(g_xi)    > 0L) sqrt(mean(g_xi^2))    else scale_mu

  # Protect against degenerate cases
  if (!is.finite(scale_mu) || scale_mu == 0) scale_mu <- 1
  ratio_sigma <- if (is.finite(scale_sigma) && scale_sigma > 0)
    scale_sigma / scale_mu else 1
  ratio_xi    <- if (is.finite(scale_xi)    && scale_xi    > 0)
    scale_xi    / scale_mu else 1

  if (verbose)
    message(sprintf(
      "  Block lambda ratios -> sigma: %.3f | xi: %.3f",
      ratio_sigma, ratio_xi
    ))

  .bic_at <- function(lam_base, init_par) {
    lam <- c(mu    = lam_base,
             sigma = lam_base * ratio_sigma,
             xi    = lam_base * ratio_xi)
    res <- .fit_at_lambda(dm, y, threshold,
                          lambda         = lam,
                          alpha          = alpha,
                          penalize_shape = penalize_shape,
                          init           = init_par,
                          maxit          = maxit,
                          calc_hessian   = FALSE,
                          obs_per_year   = obs_per_year)
    if (!res$converged || !is.finite(res$nllh_raw))
      return(list(bic = NA_real_, par = init_par, lam = lam))

    par_hat <- res$par
    # Hard threshold near-zero coefficients for BIC counting
    for (j in pen_idx_all) if (abs(par_hat[j]) < 1e-2) par_hat[j] <- 0
    k_active <- sum(par_hat[pen_idx_all] != 0) +
      (length(par_hat) - length(pen_idx_all))
    bic <- 2 * res$nllh_raw + k_active * log(n_exc)
    list(bic = bic, par = par_hat, lam = lam)
  }

  # Phase 1: coarse grid
  grid_coarse <- exp(seq(log(2000), log(0.05), length.out = 40L))
  bic_coarse  <- rep(NA_real_, length(grid_coarse))
  par_coarse  <- vector("list", length(grid_coarse))
  cur_init    <- init

  if (verbose) message("  Phase 1: coarse grid (40 points)...")

  for (i in seq_along(grid_coarse)) {
    out <- .bic_at(grid_coarse[i], cur_init)
    bic_coarse[i]  <- out$bic
    par_coarse[[i]] <- out$par
    if (!is.na(out$bic)) cur_init <- out$par

    if (verbose)
      message(sprintf("    [%2d/40] lambda=%.4f  BIC=%.2f",
                      i, grid_coarse[i],
                      if (is.na(out$bic)) NA else round(out$bic, 2)))

    # Early stopping: valley found and BIC rising consistently
    if (i >= 5L) {
      recent <- bic_coarse[max(1L, i - 3L):i]
      recent <- recent[!is.na(recent)]
      if (length(recent) >= 3L &&
          all(diff(recent) > 0) &&
          !is.na(out$bic) &&
          out$bic > min(bic_coarse[1L:(i - 1L)], na.rm = TRUE) + 50) {
        if (verbose) message("  -> Valley found, stopping coarse search early.")
        break
      }
    }
  }

  valid_c <- which(!is.na(bic_coarse))
  if (length(valid_c) == 0L)
    stop(".select_lambda_bic: all coarse grid fits failed.")

  best_c   <- valid_c[which.min(bic_coarse[valid_c])]
  lam_best <- grid_coarse[best_c]

  # Phase 2: fine grid bracketing the coarse optimum
  idx_lo     <- max(min(valid_c), best_c - 2L)
  idx_hi     <- min(max(valid_c), best_c + 2L)
  lam_lo     <- grid_coarse[idx_hi]
  lam_hi     <- grid_coarse[idx_lo]
  grid_fine  <- exp(seq(log(lam_hi), log(lam_lo), length.out = 25L))
  bic_fine   <- rep(NA_real_, length(grid_fine))
  par_fine   <- vector("list", length(grid_fine))
  cur_init   <- par_coarse[[best_c]]

  if (verbose)
    message(sprintf(
      "  Phase 2: fine grid (25 points) around lambda=%.4f...", lam_best
    ))

  for (i in seq_along(grid_fine)) {
    out <- .bic_at(grid_fine[i], cur_init)
    bic_fine[i]  <- out$bic
    par_fine[[i]] <- out$par
    if (!is.na(out$bic)) cur_init <- out$par

    if (verbose)
      message(sprintf("    [%2d/25] lambda=%.5f  BIC=%.2f",
                      i, grid_fine[i],
                      if (is.na(out$bic)) NA else round(out$bic, 2)))
  }

  valid_f <- which(!is.na(bic_fine))

  if (length(valid_f) > 0L &&
      min(bic_fine[valid_f]) < bic_coarse[best_c]) {
    best_f   <- valid_f[which.min(bic_fine[valid_f])]
    lam_best <- grid_fine[best_f]
    if (verbose)
      message(sprintf("  Fine search improved: lambda=%.5f", lam_best))
  } else {
    if (verbose)
      message(sprintf("  Coarse optimum confirmed: lambda=%.5f", lam_best))
  }

  # Return block-calibrated lambda vector
  c(mu    = lam_best,
    sigma = lam_best * ratio_sigma,
    xi    = lam_best * ratio_xi)
}
