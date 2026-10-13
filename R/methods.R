# =============================================================================
# methods.R
# Standard S3 methods for nhpp_fit objects: summary, logLik, nobs, predict,
# residuals and plot.
# =============================================================================

#' Summarise a fitted nhpp_fit model
#'
#' Builds a structured summary of a fitted \code{nhpp_fit} object, including
#' the penalty used, the active and shrunk coefficients, information
#' criteria, convergence diagnostics, and (when available) Hessian
#' conditioning. Coefficient standard errors are deliberately not reported:
#' after penalization and data-driven selection they lack a valid
#' interpretation (see \code{\link{bootstrap_coef}} for conditional intervals).
#'
#' @param object An \code{nhpp_fit} object.
#' @param tol Numeric. Coefficients smaller than this in absolute value are
#'   treated as zero. If \code{NULL}, inherits \code{object$active_tol};
#'   legacy objects fall back to \code{1e-2}.
#' @param ... Ignored.
#'
#' @return An object of class \code{summary.nhpp_fit}, printed by
#'   \code{\link{print.summary.nhpp_fit}}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' summary(fit)
#'
#' @export
summary.nhpp_fit <- function(object, tol = NULL, ...) {
  tol  <- .fit_tol(object, tol, "summary.nhpp_fit")
  par  <- object$par
  keep <- .active_parameter_mask(object, tol)

  par_oper <- if (!is.null(object$par_oper)) object$par_oper else
    ifelse(keep, par, 0)
  coefficients <- data.frame(
    estimate    = unname(par),
    operational = unname(par_oper),
    active      = unname(keep),
    row.names   = names(par)
  )

  hessian_info <- NULL
  if (!is.null(object$hessian)) {
    eigs    <- eigen(object$hessian, symmetric = TRUE,
                     only.values = TRUE)$values
    min_eig <- min(eigs)
    cond    <- max(eigs) / max(min_eig, 1e-300)
    status  <- if (any(eigs < 0)) "not_pd" else if (cond > 1e8) "ill" else "ok"
    hessian_info <- list(min_eigen = min_eig, condition = cond, status = status)
  }

  structure(
    list(
      threshold      = object$threshold,
      penalty        = object$penalty,
      lambda         = object$lambda,
      alpha          = object$alpha,
      weights_source = object$weights_source,
      penalize_shape = object$penalize_shape,
      obs_per_year   = object$obs_per_year,
      converged      = object$converged,
      kkt            = object$kkt,
      nllh_raw       = object$nllh_raw,
      nllh_pen       = object$nllh_pen,
      n_exc          = object$n_exc,
      bic            = if (is.null(object$n_exc)) NA_real_ else
        bic_nhpp(object, tol = tol),
      tol            = tol,
      coefficients   = coefficients,
      hessian        = hessian_info
    ),
    class = "summary.nhpp_fit"
  )
}

#' Print a summary of a fitted nhpp_fit model
#'
#' @param x A \code{summary.nhpp_fit} object.
#' @param digits Number of significant digits for coefficients.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)))
#' print(summary(fit_nhpp(df, threshold = 4, penalty = "none",
#'                        verbose = FALSE)))
#'
#' @export
print.summary.nhpp_fit <- function(x, digits = 5L, ...) {
  cat("-- nhpp_fit summary --------------------------------------\n")
  cat(sprintf("  Threshold      : %.4g\n",  x$threshold))
  cat(sprintf("  Penalty        : %s\n",    x$penalty))
  cat(sprintf("  Lambda         : %s\n",    .format_lambda(x$lambda)))
  if (!is.null(x$weights_source) && !identical(x$weights_source, "none"))
    cat(sprintf("  Weights        : %s\n",  x$weights_source))
  cat(sprintf("  Alpha          : %s\n",    .format_lambda(x$alpha)))
  cat(sprintf("  Penalize shape : %s\n",    x$penalize_shape))
  cat(sprintf("  obs/year       : %.2f\n",  x$obs_per_year))
  cat(sprintf("  Converged      : %s\n",    x$converged))
  if (!is.null(x$kkt) && is.finite(x$kkt))
    cat(sprintf("  KKT violation  : %.3g (per exceedance)\n", x$kkt))
  cat(sprintf("  nllh (raw)     : %.4f\n",  x$nllh_raw))
  cat(sprintf("  nllh (pen)     : %.4f\n",  x$nllh_pen))
  if (!is.null(x$n_exc))
    cat(sprintf("  Exceedances    : %d | BIC: %.3f\n", as.integer(x$n_exc), x$bic))

  co <- x$coefficients
  active <- co[co$active, , drop = FALSE]
  cat(sprintf("\n  Active coefficients (%d of %d):\n", nrow(active), nrow(co)))
  if (nrow(active) > 0L)
    print(data.frame(Estimate = round(active$estimate, digits),
                     row.names = rownames(active)))
  inactive <- rownames(co)[!co$active]
  if (length(inactive) > 0L)
    cat(sprintf("\n  Shrunk to zero : %s\n", paste(inactive, collapse = ", ")))

  if (!is.null(x$hessian)) {
    cat(sprintf("\n  Hessian: min eigenvalue = %.3g | condition = %.2e\n",
                x$hessian$min_eigen, x$hessian$condition))
    msg <- switch(x$hessian$status,
                  not_pd = "  [!] Hessian not positive definite - SEs unreliable.\n",
                  ill    = "  [!] Hessian ill-conditioned - parameters may not be identifiable.\n",
                  ok     = "  [ok] Hessian positive definite.\n")
    cat(msg)
  }
  invisible(x)
}


#' Log-likelihood, number of observations, and prediction for nhpp_fit
#'
#' \code{logLik()} returns the unpenalized log-likelihood at the operational
#' estimator, with degrees of freedom equal to the number of active
#' parameters and \code{nobs} equal to the number of exceedances, so that
#' \code{AIC()} and \code{BIC()} work; \code{BIC(fit)} equals
#' \code{\link{bic_nhpp}(fit)}. \code{predict()} is a thin wrapper around
#' \code{\link{predict_params}}.
#'
#' @param object An \code{nhpp_fit} object.
#' @param newdata Optional data frame of covariates (see
#'   \code{\link{predict_params}}).
#' @param operational Logical; see \code{\link{predict_params}}.
#' @param ... Ignored.
#'
#' @return \code{logLik()}: an object of class \code{logLik}.
#'   \code{nobs()}: the number of exceedances. \code{predict()}: a list with
#'   components \code{mu}, \code{sigma}, \code{xi}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' logLik(fit)
#' c(AIC = AIC(fit), BIC = BIC(fit), nobs = nobs(fit))
#' str(predict(fit, newdata = data.frame(x = c(-1, 0, 1))))
#'
#' @name nhpp_fit-methods
NULL

#' @rdname nhpp_fit-methods
#' @export
logLik.nhpp_fit <- function(object, ...) {
  nllh <- if (!is.null(object$nllh_oper)) object$nllh_oper else object$nllh_raw
  tol  <- .fit_tol(object, NULL, "logLik.nhpp_fit")
  structure(-nllh,
            df    = sum(.active_parameter_mask(object, tol)),
            nobs  = object$n_exc,
            class = "logLik")
}

#' @rdname nhpp_fit-methods
#' @importFrom stats nobs
#' @export
nobs.nhpp_fit <- function(object, ...) object$n_exc

#' @rdname nhpp_fit-methods
#' @export
predict.nhpp_fit <- function(object, newdata = NULL, operational = TRUE, ...) {
  predict_params(object, newdata = newdata, operational = operational)
}


#' Transformed residuals of a fitted nhpp_fit model
#'
#' Two residual types are available, both standard exponential under a
#' correctly specified model:
#' \describe{
#'   \item{\code{"interarrival"}}{Time-rescaled inter-exceedance increments
#'     \eqn{Z_k = n_y^{-1}\sum_{j: t_{k-1} \le j < t_k} \hat\tau_u(t_j)},
#'     \eqn{k = 2, \ldots, m}, where \eqn{\hat\tau_u(t)} is the fitted
#'     exceedance rate at the threshold (Ogata, 1988; Smith and Shively, 1995).}
#'   \item{\code{"magnitude"}}{Exponential transforms of the threshold
#'     excesses under the implied generalized Pareto law,
#'     \eqn{\xi^{-1}\log\{1 + \xi (y - u)/\beta(t)\}} with
#'     \eqn{\beta(t) = \sigma(t) + \xi(t)\{u - \mu(t)\}}.}
#' }
#'
#' @param object An \code{nhpp_fit} object fitted with \code{margEVT} 0.4.0
#'   or later (the observed series is stored in the object).
#' @param type Character, \code{"interarrival"} or \code{"magnitude"}.
#' @param operational Logical. Use the operational estimator (default).
#' @param ... Ignored.
#'
#' @return Numeric vector of residuals.
#'
#' @references
#' Ogata, Y. (1988). Statistical models for earthquake occurrences and
#' residual analysis for point processes. \emph{Journal of the American
#' Statistical Association}, 83(401), 9--27.
#'
#' Smith, R. L. and Shively, T. S. (1995). Point process approach to modeling
#' trends in tropospheric ozone based on exceedances of a high threshold.
#' \emph{Atmospheric Environment}, 29(23), 3489--3499.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' z <- residuals(fit)
#' c(mean = mean(z), var = var(z))
#'
#' @export
residuals.nhpp_fit <- function(object, type = c("interarrival", "magnitude"),
                               operational = TRUE, ...) {
  type <- match.arg(type)
  if (is.null(object$y))
    stop("residuals.nhpp_fit: the fitted object does not store the observed ",
         "series; refit with margEVT >= 0.4.0.")
  y   <- object$y
  u   <- object$threshold
  prm <- predict_params(object, operational = operational)
  exc <- which(y > u)

  if (type == "interarrival") {
    tau <- .tail_measure_at_level(u, prm$mu, prm$sigma, prm$xi)
    if (any(!is.finite(tau)))
      stop("residuals.nhpp_fit: the threshold lies outside the fitted support.")
    if (length(exc) < 2L) return(numeric(0L))
    cum <- c(0, cumsum(tau)) / object$obs_per_year
    # sum over t_{k-1} <= j < t_k  =  cum[t_k] - cum[t_{k-1}] (shifted index)
    return(cum[exc[-1L]] - cum[exc[-length(exc)]])
  }

  mu <- prm$mu[exc]; sig <- prm$sigma[exc]; xi <- prm$xi[exc]
  beta <- sig + xi * (u - mu)
  w <- (y[exc] - u) / beta
  ifelse(abs(xi) < 1e-6, w, log1p(xi * w) / xi)
}


#' Plot diagnostics for a fitted nhpp_fit object
#'
#' @param x An \code{nhpp_fit} object.
#' @param type Character. \code{"intensity"} (default): cumulative fitted
#'   exceedance intensity \eqn{\hat{\Lambda}(0,t)}; \code{"fitted"}:
#'   time-varying parameter profiles; \code{"path"}: BIC and coefficient
#'   paths along the regularization path (fits with \code{lambda = "bic"});
#'   \code{"qq"}: exponential Q-Q plot of transformed residuals with exact
#'   pointwise 95\% bands from the Beta distribution of exponential order
#'   statistics.
#' @param residuals Residual type for \code{type = "qq"}; see
#'   \code{\link{residuals.nhpp_fit}}.
#' @param ... Additional graphical parameters passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{NULL}, called for side effects.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' plot(fit, type = "intensity")
#' plot(fit, type = "qq")
#'
#' @export
plot.nhpp_fit <- function(x, type = c("intensity", "fitted", "path", "qq"),
                          residuals = c("interarrival", "magnitude"), ...) {
  type <- match.arg(type)

  fitted <- if (is.null(x$fitted_oper)) x$fitted else x$fitted_oper
  mu_t  <- fitted$mu
  sig_t <- fitted$sigma
  xi_t  <- fitted$xi
  u     <- x$threshold
  n_y   <- x$obs_per_year

  if (type == "intensity") {
    z_u    <- 1 + xi_t * (u - mu_t) / sig_t
    rate_t <- ifelse(
      abs(xi_t) < 1e-6,
      exp(pmax(-500, -(u - mu_t) / sig_t)),
      ifelse(z_u <= 0, 0, exp((-1 / xi_t) * log(pmax(z_u, 1e-300))))
    )
    graphics::plot(
      cumsum(rate_t) / n_y, type = "l", col = "#2C3E50", lwd = 1.5,
      main = "Integrated Intensity Measure Over Time",
      xlab = "Observation Index (t)",
      ylab = expression(hat(Lambda)(0, t)), ...
    )
    graphics::grid()

  } else if (type == "fitted") {
    var_mu  <- length(unique(round(mu_t,  6L))) > 1L
    var_sig <- length(unique(round(sig_t, 6L))) > 1L
    var_xi  <- length(unique(round(xi_t,  6L))) > 1L

    to_plot <- c(if (var_mu) "mu", if (var_sig) "sigma", if (var_xi) "xi")
    if (length(to_plot) == 0L) to_plot <- c("mu", "sigma", "xi")

    oldpar <- graphics::par(mfrow = c(length(to_plot), 1L),
                            mar = c(3.5, 4, 2, 1))
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

  } else if (type == "path") {
    if (is.null(x$path))
      stop("plot.nhpp_fit: no regularization path stored; fit with lambda = \"bic\".")
    tab <- x$path$table
    ord <- order(tab$lambda)
    tab <- tab[ord, , drop = FALSE]
    cf  <- x$path$coef[ord, , drop = FALSE]
    pen <- .penalized_parameter_indices(x)
    oldpar <- graphics::par(mfrow = c(2L, 1L), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(oldpar))
    ok <- !is.na(tab$bic)
    graphics::plot(log(tab$lambda[ok]), tab$bic[ok], type = "b", pch = 19,
                   cex = 0.6, col = "#2C3E50", xlab = expression(log(lambda)),
                   ylab = "BIC", main = "BIC along the regularization path", ...)
    graphics::abline(v = log(tab$lambda[tab$selected]), lty = 2, col = "#E74C3C")
    graphics::matplot(log(tab$lambda[ok]), cf[ok, pen, drop = FALSE],
                      type = "l", lty = 1, xlab = expression(log(lambda)),
                      ylab = "Coefficient", main = "Penalized coefficient paths")
    graphics::abline(h = 0, col = "grey60")
    graphics::abline(v = log(tab$lambda[tab$selected]), lty = 2, col = "#E74C3C")

  } else {
    residuals <- match.arg(residuals)
    z <- sort(residuals.nhpp_fit(x, type = residuals))
    m <- length(z)
    if (m < 2L) stop("plot.nhpp_fit: too few residuals for a Q-Q plot.")
    i <- seq_len(m)
    theo <- stats::qexp(i / (m + 1))
    lo <- stats::qexp(stats::qbeta(0.025, i, m - i + 1))
    hi <- stats::qexp(stats::qbeta(0.975, i, m - i + 1))
    graphics::plot(theo, z, pch = 19, cex = 0.5, col = "#2C3E50",
                   ylim = range(c(z, lo, hi)),
                   xlab = "Exp(1) quantiles", ylab = "Ordered residuals",
                   main = sprintf("Exponential Q-Q plot (%s residuals)", residuals),
                   ...)
    graphics::lines(theo, lo, lty = 2, col = "#E74C3C")
    graphics::lines(theo, hi, lty = 2, col = "#E74C3C")
    graphics::abline(0, 1, col = "grey40")
  }

  invisible(NULL)
}
