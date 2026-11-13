# =============================================================================
# print.R
# print and plot methods for the result objects other than nhpp_fit:
#   nhpp_horizon (horizon_risk), nhpp_boot (bootstrap_nhpp),
#   nhpp_backtest (backtest) and nhpp_var_generator (fit_var_generator).
# The data-frame results keep class "data.frame" as their second class, so
# indexing and column access behave as before.
# =============================================================================

#' Print and plot results of horizon_risk()
#'
#' @param x An \code{nhpp_horizon} object returned by
#'   \code{\link{horizon_risk}}.
#' @param digits Number of significant digits. Default \code{4L}.
#' @param ... Further arguments, passed to \code{\link[graphics]{matplot}} by
#'   \code{plot()} and ignored by \code{print()}.
#'
#' @return \code{print()} returns \code{x} invisibly. \code{plot()} draws, for
#'   each level, the probability that the first exceedance occurs by year
#'   \eqn{l} of the horizon, and returns \code{NULL} invisibly.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500), year = rep(1:30, each = 50))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' future <- lapply(seq(0, 2, length.out = 20),
#'                  function(v) data.frame(x = rep(v, 50)))
#' hr <- horizon_risk(fit, z = c(25, 30), trajectories = future)
#' hr
#' plot(hr)
#'
#' @name nhpp_horizon-methods
NULL

#' @rdname nhpp_horizon-methods
#' @export
print.nhpp_horizon <- function(x, digits = 4L, ...) {
  L <- attr(x, "horizon")
  S <- attr(x, "n_trajectories")
  if (!is.null(L))
    cat(sprintf("Exceedance risk over a horizon of %d year%s%s\n", L,
                if (L == 1L) "" else "s",
                if (is.null(S)) "" else
                  sprintf(" (%d trajector%s)", S, if (S == 1L) "y" else "ies")))
  print(.strip_result(x, "nhpp_horizon"), digits = digits, row.names = FALSE)
  invisible(x)
}

#' @rdname nhpp_horizon-methods
#' @export
plot.nhpp_horizon <- function(x, ...) {
  cum <- attr(x, "cumulative")
  cols <- if (is.null(cum)) NA else
    match(x$z, suppressWarnings(as.numeric(colnames(cum))))
  if (is.null(cum) || length(cols) == 0L || anyNA(cols))
    stop("plot.nhpp_horizon: no cumulative probabilities stored for these ",
         "levels; plot the result of horizon_risk() or rows of it.")
  cum <- cum[, cols, drop = FALSE]
  L <- nrow(cum)
  graphics::matplot(seq_len(L), cum, type = "l", lty = 1, lwd = 1.5,
                    col = .result_palette(ncol(cum)), ylim = c(0, 1),
                    xlab = "Year of the horizon",
                    ylab = "P(first exceedance by year l)",
                    main = "Exceedance risk over the horizon", ...)
  graphics::grid()
  graphics::legend("topleft", legend = paste("z =", colnames(cum)),
                   col = .result_palette(ncol(cum)), lty = 1, lwd = 1.5,
                   bty = "n")
  invisible(NULL)
}


#' Print results of bootstrap_nhpp()
#'
#' @param x An \code{nhpp_boot} object returned by
#'   \code{\link{bootstrap_nhpp}}.
#' @param digits Number of significant digits. Default \code{5L}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500), year = rep(1:30, each = 50))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' \donttest{
#' bootstrap_nhpp(fit, df, TRs = c(10, 50), R = 30, verbose = FALSE)
#' }
#'
#' @export
print.nhpp_boot <- function(x, digits = 5L, ...) {
  info <- attr(x, "info")
  cat("Conditional parametric bootstrap")
  if (!is.null(info))
    cat(sprintf(" | R = %d | %.0f%% percentile intervals | occurrence: %s",
                info$R, 100 * info$level, info$occurrence))
  cat("\n")
  if (!is.null(x$rl)) {
    cat(sprintf("\nReturn levels (approach %s):\n",
                if (is.null(info)) "?" else info$approach))
    print(x$rl, digits = digits, row.names = FALSE)
  }
  cat("\nCoefficients:\n")
  print(x$coef, digits = digits, row.names = FALSE)
  invisible(x)
}


#' Print results of backtest()
#'
#' Prints the nominal binomial comparisons of the exceedance rates and the
#' exploratory tail-PIT Kolmogorov-Smirnov comparison. The per-year results
#' are in \code{x$results}.
#'
#' @param x An \code{nhpp_backtest} object returned by
#'   \code{\link{backtest}}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(30 * 50), year = rep(1991:2020, each = 50))
#' df$y <- 5 + df$x + rexp(nrow(df), 1 / 2)
#' fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' \donttest{
#' bt <- backtest(fit, df, varname = "y", TRs = c(2, 5), n_obs = 50,
#'                window_years = 5, min_train_years = 15, n_boot = 50,
#'                verbose = FALSE)
#' bt
#' }
#'
#' @export
print.nhpp_backtest <- function(x, ...) {
  n_years <- if (is.null(x$results)) 0L else nrow(x$results)
  cat(sprintf("Walk-forward backtest: %d validation year%s\n", n_years,
              if (n_years == 1L) "" else "s"))
  if (!is.null(x$binom_tests)) {
    cat("\nExceedance rates (binomial comparison with 1/T):\n")
    b <- x$binom_tests
    b$calibrated <- NULL
    print(b, row.names = FALSE)
  }
  if (!is.null(x$ks_test)) {
    cat("\nTail-conditional PIT, Kolmogorov-Smirnov (exploratory):\n")
    k <- x$ks_test
    k$calibrated <- NULL
    print(k, row.names = FALSE)
  }
  invisible(x)
}


#' Print a fitted VAR covariate generator
#'
#' @param x An \code{nhpp_var_generator} object returned by
#'   \code{\link{fit_var_generator}}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' n <- 3650
#' t <- rep(1:365, 10)
#' df <- cbind(data.frame(x = as.numeric(stats::arima.sim(list(ar = 0.6), n))),
#'             seasonal_harmonics(t))
#' df$y <- 10 + df$x + rexp(n, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 365, verbose = FALSE)
#' fit_var_generator(fit, df)
#'
#' @export
print.nhpp_var_generator <- function(x, ...) {
  cat(sprintf("VAR(%d) covariate generator\n", x$p_opt))
  cat(sprintf("  Covariates      : %s\n", paste(x$vars, collapse = ", ")))
  det <- setdiff(x$vars, x$stochastic_vars)
  if (length(det))
    cat(sprintf("  Deterministic   : %s (seasonal cycle only)\n",
                paste(det, collapse = ", ")))
  cat(sprintf("  Seasonal model  : %d harmonic%s, period %s\n",
              x$n_harmonics, if (x$n_harmonics == 1L) "" else "s",
              format(x$period)))
  if (x$K > 0L)
    cat(sprintf("  Spectral radius : %.4f (stable)\n", x$spectral_radius))
  if (is.function(x$derive))
    cat(sprintf("  Derived         : model covariates built by `derive`%s\n",
                if (isTRUE(x$derive_lead > 0L))
                  sprintf(" (lead of %d observations)", x$derive_lead) else ""))
  invisible(x)
}


# --- helpers -----------------------------------------------------------------

# Drops the result class and the bulky attributes before printing.
.strip_result <- function(x, cls) {
  class(x) <- setdiff(class(x), cls)
  for (a in c("cumulative", "horizon", "n_trajectories")) attr(x, a) <- NULL
  x
}

.result_palette <- function(n) {
  grDevices::hcl.colors(max(n, 2L), "Dark 3")[seq_len(n)]
}
