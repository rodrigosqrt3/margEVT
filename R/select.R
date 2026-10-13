# =============================================================================
# select.R
# BIC-based lambda selection for penalized NHPP models.
#
# .select_lambda_bic() is internal - called by fit_nhpp() when lambda = "bic".
#
# Algorithm:
#   Default: Phase 1 coarse log-spaced grid (early stopping when a BIC valley
#   is found), then Phase 2 fine grid bracketing the coarse optimum.
#   With a user-supplied `lambda_grid`, a single pass over that grid is made.
#   Each path solution is warm-started from the previous smooth solution; the
#   BIC is evaluated at the operational (activity-thresholded) estimator.
#   The selected lambda is returned together with its path solution, so the
#   final model is exactly the one scored by the criterion.
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
#' @param active_tol Positive numeric tolerance for BIC complexity counting.
#' @param lambda_scaling Character, \code{"common"} or \code{"gradient"}.
#' @param maxit Integer.
#' @param verbose Logical.
#' @param penalty_weights Optional coordinate weights (see \code{fit_nhpp}).
#' @param lambda_grid Optional user grid of base penalties.
#' @param grad_tol KKT tolerance (see \code{fit_nhpp}).
#'
#' @return A list with \code{lambda} (named vector \code{c(mu=, sigma=, xi=)}),
#'   \code{par} (smooth path solution at the selected lambda), and
#'   \code{path} (a list with a data frame \code{table} and a coefficient
#'   matrix \code{coef}).
#'
#' @noRd
.select_lambda_bic <- function(dm, y, threshold, alpha,
                               penalize_shape, init,
                               obs_per_year    = 365.25,
                               active_tol      = 1e-2,
                               lambda_scaling  = c("common", "gradient"),
                               maxit           = 5000L,
                               verbose         = TRUE,
                               penalty_weights = NULL,
                               lambda_grid     = NULL,
                               grad_tol        = 5e-4) {

  lambda_scaling <- match.arg(lambda_scaling)

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  n_exc <- sum(y > threshold, na.rm = TRUE)

  pen_idx_all <- .penalized_index(dm, penalize_shape)

  ratio_sigma <- ratio_xi <- 1
  if (lambda_scaling == "gradient") {
    # Scale each block's lambda by its score magnitude at the initial model.
    # This is a numerical calibration, not an adaptive-LASSO weight.
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

    if (!is.finite(scale_mu) || scale_mu == 0) scale_mu <- 1
    ratio_sigma <- if (is.finite(scale_sigma) && scale_sigma > 0)
      scale_sigma / scale_mu else 1
    ratio_xi <- if (is.finite(scale_xi) && scale_xi > 0)
      scale_xi / scale_mu else 1
  }

  if (verbose)
    message(sprintf(
      "  Lambda scaling: %s | ratios -> sigma: %.3f | xi: %.3f",
      lambda_scaling, ratio_sigma, ratio_xi
    ))

  block_lambda <- function(lam_base)
    c(mu = lam_base, sigma = lam_base * ratio_sigma, xi = lam_base * ratio_xi)

  .bic_at <- function(lam_base, init_par) {
    lam <- block_lambda(lam_base)
    res <- .fit_at_lambda(dm, y, threshold,
                          lambda          = lam,
                          alpha           = alpha,
                          penalize_shape  = penalize_shape,
                          init            = init_par,
                          maxit           = maxit,
                          calc_hessian    = FALSE,
                          obs_per_year    = obs_per_year,
                          penalty_weights = penalty_weights,
                          grad_tol        = grad_tol,
                          active_tol      = active_tol)
    if (!res$converged || !is.finite(res$nllh_raw))
      return(list(bic = NA_real_, par_smooth = init_par,
                  par_oper = init_par, k_active = NA_integer_,
                  converged = FALSE))

    par_smooth <- res$par
    par_oper   <- par_smooth
    # The operational estimator defines both the BIC likelihood and its
    # effective parameter count. The smooth optimizer solution remains the
    # warm start for the next lambda on the path.
    for (j in pen_idx_all)
      if (abs(par_oper[j]) < active_tol) par_oper[j] <- 0
    nllh_oper <- pp_nllh(par_oper, dm, y, threshold,
                         lambda = 0, alpha = alpha,
                         pen_xi = penalize_shape,
                         obs_per_year = obs_per_year)
    if (!is.finite(nllh_oper) || nllh_oper >= 1e9)
      return(list(bic = NA_real_, par_smooth = init_par,
                  par_oper = init_par, k_active = NA_integer_,
                  converged = FALSE))
    k_active <- sum(par_oper[pen_idx_all] != 0) +
      (length(par_oper) - length(pen_idx_all))
    list(bic = 2 * nllh_oper + k_active * log(n_exc),
         par_smooth = par_smooth, par_oper = par_oper,
         k_active = as.integer(k_active), converged = TRUE)
  }

  rows <- list()
  record <- function(lam_base, out, phase) {
    rows[[length(rows) + 1L]] <<- list(
      lambda = lam_base, phase = phase, bic = out$bic,
      k_active = out$k_active, converged = out$converged,
      par_smooth = out$par_smooth, par_oper = out$par_oper)
  }

  run_grid <- function(grid, start, phase, early_stop) {
    bic <- rep(NA_real_, length(grid))
    pars <- vector("list", length(grid))
    cur <- start
    for (i in seq_along(grid)) {
      out <- .bic_at(grid[i], cur)
      bic[i] <- out$bic
      pars[[i]] <- out$par_smooth
      record(grid[i], out, phase)
      if (!is.na(out$bic)) cur <- out$par_smooth
      if (verbose)
        message(sprintf("    [%2d/%d] lambda=%.5f  BIC=%.2f",
                        i, length(grid), grid[i],
                        if (is.na(out$bic)) NA else round(out$bic, 2)))
      # Early stopping: valley found and BIC rising consistently
      if (early_stop && i >= 5L) {
        recent <- bic[max(1L, i - 3L):i]
        recent <- recent[!is.na(recent)]
        if (length(recent) >= 3L &&
            all(diff(recent) > 0) &&
            !is.na(out$bic) &&
            out$bic > min(bic[1L:(i - 1L)], na.rm = TRUE) + 50) {
          if (verbose) message("  -> Valley found, stopping coarse search early.")
          break
        }
      }
    }
    list(bic = bic, pars = pars)
  }

  if (!is.null(lambda_grid)) {
    grid <- sort(unique(lambda_grid), decreasing = TRUE)
    if (verbose) message(sprintf("  User grid (%d points)...", length(grid)))
    g <- run_grid(grid, init, "user", early_stop = FALSE)
    valid <- which(!is.na(g$bic))
    if (length(valid) == 0L)
      stop(".select_lambda_bic: all grid fits failed.")
    best <- valid[which.min(g$bic[valid])]
    lam_best <- grid[best]
    par_best <- g$pars[[best]]
  } else {
    # Phase 1: coarse grid
    grid_coarse <- exp(seq(log(2000), log(0.05), length.out = 40L))
    if (verbose) message("  Phase 1: coarse grid (40 points)...")
    gc <- run_grid(grid_coarse, init, "coarse", early_stop = TRUE)
    valid_c <- which(!is.na(gc$bic))
    if (length(valid_c) == 0L)
      stop(".select_lambda_bic: all coarse grid fits failed.")
    best_c   <- valid_c[which.min(gc$bic[valid_c])]
    lam_best <- grid_coarse[best_c]
    par_best <- gc$pars[[best_c]]

    # Phase 2: fine grid bracketing the coarse optimum
    idx_lo    <- max(min(valid_c), best_c - 2L)
    idx_hi    <- min(max(valid_c), best_c + 2L)
    grid_fine <- exp(seq(log(grid_coarse[idx_lo]), log(grid_coarse[idx_hi]),
                         length.out = 25L))
    if (verbose)
      message(sprintf(
        "  Phase 2: fine grid (25 points) around lambda=%.4f...", lam_best
      ))
    gf <- run_grid(grid_fine, par_best, "fine", early_stop = FALSE)
    valid_f <- which(!is.na(gf$bic))
    if (length(valid_f) > 0L &&
        min(gf$bic[valid_f]) < gc$bic[best_c]) {
      best_f   <- valid_f[which.min(gf$bic[valid_f])]
      lam_best <- grid_fine[best_f]
      par_best <- gf$pars[[best_f]]
      if (verbose)
        message(sprintf("  Fine search improved: lambda=%.5f", lam_best))
    } else if (verbose) {
      message(sprintf("  Coarse optimum confirmed: lambda=%.5f", lam_best))
    }
  }

  table <- data.frame(
    lambda    = vapply(rows, `[[`, numeric(1L), "lambda"),
    phase     = vapply(rows, `[[`, character(1L), "phase"),
    bic       = vapply(rows, `[[`, numeric(1L), "bic"),
    k_active  = vapply(rows, function(r) as.integer(r$k_active), integer(1L)),
    converged = vapply(rows, `[[`, logical(1L), "converged"),
    stringsAsFactors = FALSE
  )
  table$selected <- FALSE
  cand <- which(table$lambda == lam_best & !is.na(table$bic))
  table$selected[cand[which.min(table$bic[cand])]] <- TRUE
  coef_path <- do.call(rbind, lapply(rows, `[[`, "par_oper"))
  colnames(coef_path) <- names(init)

  list(lambda = block_lambda(lam_best),
       par    = par_best,
       path   = list(table = table, coef = coef_path))
}
