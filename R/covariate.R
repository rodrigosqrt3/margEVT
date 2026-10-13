# =============================================================================
# covariate.R
# Build covariate data frames for marginalization.
#
# The central function is build_cov_annual(): given a fitted model and a
# named list of covariate values, it returns a data frame with one row per
# time point (typically 365 rows for one year of daily data) containing
# all columns the model needs.
#
# Seasonal harmonics are recognised by name: cos<k>, sin<k>, and the legacy
# spelling sen<k> (k = 1, 2, ...). They are evaluated at positions
# t = 1, ..., n_obs within the year, using the same convention as
# seasonal_harmonics().
# =============================================================================

# Value of one named harmonic column at time positions `t`.
.harmonic_column <- function(name, t, period) {
  k <- as.integer(sub("^(cos|sin|sen)", "", name))
  if (startsWith(name, "cos")) cos(2 * pi * k * t / period) else
    sin(2 * pi * k * t / period)
}

#' Seasonal Fourier harmonics
#'
#' Builds the seasonal harmonic columns used by \code{margEVT}. The
#' \eqn{k}-th pair is \eqn{\cos(2\pi k t / P)} and \eqn{\sin(2\pi k t / P)},
#' where \eqn{t} is the position of the observation within the year and
#' \eqn{P} is the seasonal period. Using this function to build the fitting
#' data guarantees that \code{\link{build_cov_annual}} reconstructs the same
#' seasonal cycle during marginalization, where positions run from 1 to
#' \code{n_obs} within each simulated or resampled year.
#'
#' @param t Numeric vector of within-year positions (e.g. day of year).
#' @param period Positive seasonal period. Default \code{365.25}.
#' @param n_harmonics Positive integer number of harmonic pairs. Default 2.
#' @param sine Prefix of the sine columns: \code{"sen"} (default, the
#'   historical margEVT naming) or \code{"sin"}.
#'
#' @return A data frame with columns \code{cos1}, \code{sen1} (or
#'   \code{sin1}), ..., one row per element of \code{t}.
#'
#' @examples
#' h <- seasonal_harmonics(1:365, period = 365.25)
#' head(h)
#' names(seasonal_harmonics(1:52, period = 52, n_harmonics = 3, sine = "sin"))
#'
#' @export
seasonal_harmonics <- function(t, period = 365.25, n_harmonics = 2L,
                               sine = c("sen", "sin")) {
  sine <- match.arg(sine)
  if (!is.numeric(t) || any(!is.finite(t)))
    stop("seasonal_harmonics: `t` must be a finite numeric vector.")
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) ||
      period <= 0)
    stop("seasonal_harmonics: `period` must be a single positive number.")
  if (!is.numeric(n_harmonics) || length(n_harmonics) != 1L ||
      !is.finite(n_harmonics) || n_harmonics < 1)
    stop("seasonal_harmonics: `n_harmonics` must be a positive integer.")
  out <- list()
  for (k in seq_len(as.integer(n_harmonics))) {
    out[[paste0("cos", k)]] <- cos(2 * pi * k * t / period)
    out[[paste0(sine, k)]]  <- sin(2 * pi * k * t / period)
  }
  as.data.frame(out)
}

#' Build a covariate data frame for marginalization
#'
#' Constructs a data frame with \code{n_obs} rows containing all columns
#' required by a fitted \code{nhpp_fit} model. The default seasonal columns
#' (\code{cos1}, \code{sen1}, \code{cos2}, \code{sen2}) and any other
#' harmonic column required by the model (\code{cos<k>}, \code{sin<k>},
#' \code{sen<k>}) are computed automatically with \code{period} (see
#' \code{\link{seasonal_harmonics}}). Any other required column is filled from
#' \code{cov_vals}, or set to zero with a warning if absent.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param cov_vals Named list of covariate values. Each element is either
#'   a scalar (replicated to \code{n_obs} rows) or a vector of length
#'   \code{n_obs}.
#' @param n_obs Integer number of rows (time points) to generate. If
#'   \code{NULL}, uses the rounded \code{fit$obs_per_year} value.
#' @param period Numeric seasonal period in observations. If \code{NULL},
#'   inherits \code{fit$obs_per_year}; legacy objects fall back to \code{365.25}.
#' @param interactions Named list of length-2 character vectors declaring
#'   interaction columns to compute. Each element is
#'   \code{c("col_A", "col_B")} and the interaction column is their product.
#'   The name of the list element becomes the column name.
#'   Example: \code{list(AB = c("A", "B"))} adds column \code{AB = A * B}.
#'
#' @return A data frame with \code{n_obs} rows and all columns required
#'   by \code{fit}.
#'
#' @examples
#' set.seed(1)
#' n <- 3650
#' df <- cbind(data.frame(x = rnorm(n)), seasonal_harmonics(rep(1:365, 10)))
#' df$y <- 10 + df$x + 3 * df$cos1 + rexp(n, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 20, loc_vars = c("x", "cos1", "sen1"),
#'                 free_vars = c("cos1", "sen1"), penalty = "none",
#'                 verbose = FALSE)
#' head(build_cov_annual(fit, list(x = 0.5), n_obs = 365))
#'
#' @export
build_cov_annual <- function(fit, cov_vals = list(), n_obs = NULL,
                             period = NULL, interactions = list()) {

  if (!inherits(fit, "nhpp_fit"))
    stop("build_cov_annual: `fit` must be an nhpp_fit object.")
  if (!is.list(interactions) ||
      (length(interactions) > 0L &&
       (is.null(names(interactions)) || any(names(interactions) == "") ||
        anyDuplicated(names(interactions)))))
    stop("build_cov_annual: `interactions` must be a named list with unique names.")
  invalid_interactions <- names(interactions)[
    vapply(interactions, length, integer(1L)) != 2L
  ]
  if (length(invalid_interactions) > 0L)
    stop(sprintf(
      "build_cov_annual: interaction '%s' must name exactly 2 columns.",
      invalid_interactions[1L]
    ))

  n_obs <- if (is.null(n_obs)) {
    as.integer(round(if (is.null(fit$obs_per_year)) 365.25 else fit$obs_per_year))
  } else as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L)
    stop("build_cov_annual: `n_obs` must be a positive integer.")
  period <- if (is.null(period)) {
    if (is.null(fit$obs_per_year)) 365.25 else fit$obs_per_year
  } else period
  if (!is.numeric(period) || length(period) != 1L ||
      !is.finite(period) || period <= 0)
    stop("build_cov_annual: `period` must be a single positive numeric value.")

  dm        <- fit$dm
  need_cols <- unique(c(
    colnames(dm$X_mu),
    colnames(dm$X_sigma),
    colnames(dm$X_xi)
  ))
  need_cols <- need_cols[need_cols != "(Intercept)"]

  t_grid <- seq_len(n_obs)
  harmonic_cols <- unique(c("cos1", "sen1", "cos2", "sen2",
                            need_cols[.is_harmonic(need_cols)]))
  df <- as.data.frame(stats::setNames(
    lapply(harmonic_cols, .harmonic_column, t = t_grid, period = period),
    harmonic_cols))

  interaction_output_cols <- names(interactions)

  add_covariate <- function(v) {
    val <- cov_vals[[v]]
    if (!is.atomic(val) || !is.numeric(val) ||
        !length(val) %in% c(1L, n_obs) || any(!is.finite(val)))
      stop(sprintf(
        paste0(
          "build_cov_annual: covariate '%s' must be a finite numeric ",
          "scalar or a vector of length n_obs (%d)."
        ),
        v, n_obs
      ))
    df[[v]] <<- if (length(val) == n_obs) val else rep(val, n_obs)
  }

  interaction_inputs <- unique(unlist(interactions, use.names = FALSE))
  for (v in interaction_inputs) {
    if (v %in% names(df)) next
    if (!v %in% names(cov_vals))
      stop(sprintf(
        "build_cov_annual: interaction input '%s' is absent from `cov_vals`.", v
      ))
    add_covariate(v)
  }

  missing_covariates <- character(0L)
  for (v in need_cols) {
    if (v %in% names(df)) next
    if (v %in% interaction_output_cols) next
    if (v %in% names(cov_vals)) {
      add_covariate(v)
    } else {
      df[[v]] <- 0
      missing_covariates <- c(missing_covariates, v)
    }
  }

  if (length(missing_covariates) > 0L)
    warning(sprintf(
      "build_cov_annual: column(s) %s not in `cov_vals`; filled with 0.",
      paste(sprintf("'%s'", missing_covariates), collapse = ", ")
    ))

  for (nm in names(interactions)) {
    cols    <- interactions[[nm]]
    df[[nm]] <- df[[cols[1L]]] * df[[cols[2L]]]
  }

  df
}


#' Extract names of non-seasonal covariates from a fitted model
#'
#' Returns the covariate names that appear in the model with a non-negligible
#' coefficient - i.e. those that actually affect the parameters and need to
#' be supplied during marginalization. Seasonal harmonic columns
#' (\code{cos<k>}, \code{sin<k>}, \code{sen<k>}) are always excluded.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param free_cols Character vector of column names to always exclude from
#'   the result (intercept and seasonality terms are always excluded).
#'   Default \code{NULL}.
#' @param tol Numeric. Coefficients smaller than this in absolute value are
#'   treated as zero. If \code{NULL}, inherits \code{fit$active_tol}; legacy
#'   objects fall back to \code{1e-2}.
#'
#' @return Character vector of active covariate names, or \code{NULL} if
#'   the model has no active non-seasonal covariates.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(2000), z = rnorm(2000))
#' df$y <- 10 + 2 * df$x + rexp(2000, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = c("x", "z"),
#'                 penalty = "lasso", lambda = 5, verbose = FALSE)
#' active_covariates(fit)
#'
#' @export
active_covariates <- function(fit, free_cols = NULL, tol = NULL) {

  if (!inherits(fit, "nhpp_fit"))
    stop("active_covariates: `fit` must be an nhpp_fit object.")
  tol <- .fit_tol(fit, tol, "active_covariates")

  dm  <- fit$dm
  par <- fit$par

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  beta <- list(
    stats::setNames(par[seq_len(p_mu)], colnames(dm$X_mu)),
    stats::setNames(par[p_mu + seq_len(p_sig)], colnames(dm$X_sigma)),
    stats::setNames(par[p_mu + p_sig + seq_len(p_xi)], colnames(dm$X_xi))
  )

  active <- character(0L)
  for (b in beta) {
    cn <- names(b)
    cn <- cn[cn != "(Intercept)" & !.is_harmonic(cn) & !cn %in% free_cols]
    active <- c(active, cn[abs(b[cn]) > tol])
  }

  active <- unique(active)
  if (length(active) == 0L) NULL else active
}
