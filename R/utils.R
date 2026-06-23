# =============================================================================
# utils.R
# Small shared utilities.
# =============================================================================

#' Summarise a fitted nhpp_fit model
#'
#' Prints a structured summary of a fitted \code{nhpp_fit} object, including
#' the penalty used, active coefficients, and basic model diagnostics.
#'
#' @param object An \code{nhpp_fit} object.
#' @param tol Numeric. Coefficients smaller than this in absolute value are
#'   treated as zero. Default \code{1e-4}.
#' @param ... Ignored.
#'
#' @export
summary.nhpp_fit <- function(object, tol = 1e-4, ...) {
  cat("-- nhpp_fit summary --------------------------------------\n")
  cat(sprintf("  Threshold      : %.4g\n",  object$threshold))
  cat(sprintf("  Penalty        : %s\n",    object$penalty))
  cat(sprintf("  Lambda (mean)  : %.5g\n",  mean(object$lambda)))
  cat(sprintf("  Alpha          : %.3g\n",  mean(object$alpha)))
  cat(sprintf("  Penalize shape : %s\n",    object$penalize_shape))
  cat(sprintf("  obs/year       : %.2f\n",  object$obs_per_year))
  cat(sprintf("  Converged      : %s\n",    object$converged))
  cat(sprintf("  nllh (raw)     : %.4f\n",  object$nllh_raw))
  cat(sprintf("  nllh (pen)     : %.4f\n",  object$nllh_pen))

  par      <- object$par
  active   <- par[abs(par) >  tol]
  inactive <- par[abs(par) <= tol]

  cat(sprintf("\n  Active coefficients (%d of %d):\n",
              length(active), length(par)))
  if (length(active) > 0L) {
    df_act <- data.frame(
      Estimate = round(active, 5L),
      row.names = names(active)
    )
    print(df_act)
  }

  if (length(inactive) > 0L)
    cat(sprintf("\n  Shrunk to zero : %s\n",
                paste(names(inactive), collapse = ", ")))

  if (!is.null(object$hessian)) {
    eigs    <- eigen(object$hessian, symmetric = TRUE,
                     only.values = TRUE)$values
    min_eig <- min(eigs)
    cond    <- max(eigs) / max(min_eig, 1e-300)
    cat(sprintf("\n  Hessian: min eigenvalue = %.3g | condition = %.2e\n",
                min_eig, cond))
    if (any(eigs < 0))
      cat("  [!] Hessian not positive definite - SEs unreliable.\n")
    else if (cond > 1e8)
      cat("  [!] Hessian ill-conditioned - parameters may not be identifiable.\n")
    else
      cat("  [ok] Hessian positive definite.\n")
  }

  invisible(object)
}


#' Compute BIC for a fitted nhpp_fit model
#'
#' @param fit An \code{nhpp_fit} object.
#' @param tol Numeric. Coefficients smaller than this are counted as zero
#'   for the active parameter count. Default \code{1e-2}.
#'
#' @return Numeric scalar. BIC value.
#'
#' @export
bic_nhpp <- function(fit, tol = 1e-2) {
  if (!inherits(fit, "nhpp_fit"))
    stop("bic_nhpp: `fit` must be an nhpp_fit object.")
  if (!is.finite(fit$nllh_raw))
    return(NA_real_)
  n_exc    <- sum(fit$dm$X_mu[, 1L])   # proxy: nrow of training data
  # Better: count exceedances implied by threshold
  # We don't store n_exc, so use active parameter count only
  k_active <- sum(abs(fit$par) > tol)
  2 * fit$nllh_raw + k_active * log(max(k_active, 2L))
}


#' Check if an object is an nhpp_fit
#'
#' @param x Any R object.
#' @return Logical.
#' @export
is_nhpp_fit <- function(x) inherits(x, "nhpp_fit")


#' Extract the number of exceedances from a fitted model
#'
#' @param fit An \code{nhpp_fit} object.
#' @param y Numeric vector of observations used to fit the model.
#'
#' @return Integer. Number of observations exceeding the threshold.
#'
#' @export
n_exceedances <- function(fit, y) {
  if (!inherits(fit, "nhpp_fit"))
    stop("n_exceedances: `fit` must be an nhpp_fit object.")
  sum(y > fit$threshold, na.rm = TRUE)
}


#' Wide-format return level table from marginalize() output
#'
#' Converts the long-format output of \code{\link{marginalize}} into a
#' wide table with one row per scenario and one column per return period.
#'
#' @param marg_result Data frame returned by \code{\link{marginalize}}.
#'
#' @return A wide data frame.
#'
#' @export
rl_table <- function(marg_result) {
  if (!is.data.frame(marg_result))
    stop("rl_table: `marg_result` must be a data frame from marginalize().")
  if (!all(c("approach", "scenario", "TR", "RL") %in% names(marg_result)))
    stop("rl_table: expected columns: approach, scenario, TR, RL.")

  keys     <- unique(marg_result[, c("approach", "scenario")])
  tr_vals  <- sort(unique(marg_result$TR))
  tr_cols  <- paste0("T", tr_vals)

  result <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    app <- keys$approach[i]
    sc  <- keys$scenario[i]
    sub <- marg_result[marg_result$approach == app &
                         marg_result$scenario == sc, ]
    rl_vals <- sub$RL[match(tr_vals, sub$TR)]
    row <- data.frame(approach = app, scenario = sc)
    for (j in seq_along(tr_cols)) row[[tr_cols[j]]] <- rl_vals[j]
    row
  }))

  result
}
