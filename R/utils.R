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
#' @return Invisibly returns the input \code{nhpp_fit} object.
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

  k_active <- sum(abs(fit$par) > tol)
  2 * fit$nllh_raw + k_active * log(fit$n_exc)
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

#' Plot diagnostics for a fitted nhpp_fit object
#'
#' Generates primary diagnostic plots for a fitted \code{nhpp_fit} model,
#' including the cumulative intensity measure over time and parameter profiles.
#'
#' @param x An \code{nhpp_fit} object.
#' @param type Character. Type of plot to generate: \code{"intensity"} (default)
#'   for the cumulative intensity measure \eqn{\hat{\Lambda}(0,t)}, or
#'   \code{"fitted"} for time-varying parameter profiles (\eqn{\mu(t)}, \eqn{\sigma(t)},
#'   and \eqn{\xi(t)}).
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{NULL}, called for side effects.
#'
#' @export
plot.nhpp_fit <- function(x, type = c("intensity", "fitted"), ...) {
  type <- match.arg(type)

  mu_t  <- x$fitted$mu
  sig_t <- x$fitted$sigma
  xi_t  <- x$fitted$xi
  u     <- x$threshold
  n_y   <- x$obs_per_year

  if (type == "intensity") {
    z_u    <- 1 + xi_t * (u - mu_t) / sig_t
    rate_t <- ifelse(
      abs(xi_t) < 1e-6,
      exp(pmax(-500, -(u - mu_t) / sig_t)),
      ifelse(z_u <= 0, 0, exp((-1 / xi_t) * log(pmax(z_u, 1e-300))))
    )

    cum_intensity <- cumsum(rate_t) / n_y

    graphics::plot(
      cum_intensity, type = "l", col = "#2C3E50", lwd = 1.5,
      main = "Integrated Intensity Measure Over Time",
      xlab = "Observation Index (t)",
      ylab = expression(hat(Lambda)(0, t)), ...
    )
    graphics::grid()

  } else {
    var_mu  <- length(unique(round(mu_t,  6L))) > 1L
    var_sig <- length(unique(round(sig_t, 6L))) > 1L
    var_xi  <- length(unique(round(xi_t,  6L))) > 1L

    to_plot <- c(if (var_mu) "mu", if (var_sig) "sigma", if (var_xi) "xi")
    if (length(to_plot) == 0L) {
      to_plot <- c("mu", "sigma", "xi")
    }

    n_panels <- length(to_plot)
    oldpar   <- graphics::par(mfrow = c(n_panels, 1L), mar = c(3.5, 4, 2, 1))
    on.exit(graphics::par(oldpar))

    if ("mu" %in% to_plot) {
      graphics::plot(mu_t, type = "l", col = "#2C3E50", lwd = 1.2,
                     main = expression("Fitted Location Parameter " * mu(t)),
                     xlab = "Time Index", ylab = expression(mu(t)), ...)
      graphics::grid()
    }

    if ("sigma" %in% to_plot) {
      graphics::plot(sig_t, type = "l", col = "#E74C3C", lwd = 1.2,
                     main = expression("Fitted Scale Parameter " * sigma(t)),
                     xlab = "Time Index", ylab = expression(sigma(t)), ...)
      graphics::grid()
    }

    if ("xi" %in% to_plot) {
      graphics::plot(xi_t, type = "l", col = "#27AE60", lwd = 1.2,
                     main = expression("Fitted Shape Parameter " * xi(t)),
                     xlab = "Time Index", ylab = expression(xi(t)), ...)
      graphics::grid()
    }
  }

  invisible(NULL)
}
