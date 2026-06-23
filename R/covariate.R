# =============================================================================
# covariate.R
# Build covariate data frames for marginalization.
#
# The central function is build_cov_annual(): given a fitted model and a
# named list of covariate values, it returns a data frame with one row per
# time point (typically 365 rows for one year of daily data) containing
# all columns the model needs.
#
# No variable names are hardcoded. Interactions must be declared explicitly
# by the user via the `interactions` argument.
# =============================================================================

#' Build a covariate data frame for marginalization
#'
#' Constructs a data frame with \code{n_obs} rows containing all columns
#' required by a fitted \code{nhpp_fit} model. Seasonality columns
#' (\code{cos1}, \code{sen1}, \code{cos2}, \code{sen2}) are computed
#' automatically from \code{n_obs}. Any other required column is filled
#' from \code{cov_vals}, or set to zero with a warning if absent.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param cov_vals Named list of covariate values. Each element is either
#'   a scalar (replicated to \code{n_obs} rows) or a vector of length
#'   \code{n_obs}.
#' @param n_obs Integer. Number of rows (time points) to generate.
#'   Typically \code{365L} for one year of daily data.
#' @param period Numeric. Period for the seasonal harmonics in days.
#'   Default \code{365.25}.
#' @param interactions Named list of length-2 character vectors declaring
#'   interaction columns to compute. Each element is
#'   \code{c("col_A", "col_B")} and the interaction column is their product.
#'   The name of the list element becomes the column name.
#'   Example: \code{list(AB = c("A", "B"))} adds column \code{AB = A * B}.
#'
#' @return A data frame with \code{n_obs} rows and all columns required
#'   by \code{fit}.
#'
#' @export
build_cov_annual <- function(fit, cov_vals = list(), n_obs = 365L,
                             period = 365.25, interactions = list()) {

  if (!inherits(fit, "nhpp_fit"))
    stop("build_cov_annual: `fit` must be an nhpp_fit object.")

  n_obs <- as.integer(n_obs)

  t_grid <- seq_len(n_obs)
  df <- data.frame(
    cos1 = cos(2 * pi * t_grid / period),
    sen1 = sin(2 * pi * t_grid / period),
    cos2 = cos(4 * pi * t_grid / period),
    sen2 = sin(4 * pi * t_grid / period)
  )

  dm        <- fit$dm
  need_cols <- unique(c(
    colnames(dm$X_mu),
    colnames(dm$X_sigma),
    colnames(dm$X_xi)
  ))
  need_cols <- need_cols[need_cols != "(Intercept)"]

  interaction_output_cols <- names(interactions)

  known_at_this_point <- c(names(df), names(cov_vals))
  for (nm in names(interactions)) {
    cols <- interactions[[nm]]
    if (length(cols) != 2L)
      stop(sprintf(
        "build_cov_annual: interaction '%s' must name exactly 2 columns.", nm
      ))
    missing_inter <- setdiff(cols, known_at_this_point)
    if (length(missing_inter) > 0L)
      stop(sprintf(
        "build_cov_annual: interaction '%s' requires column(s) not available: %s",
        nm, paste(missing_inter, collapse = ", ")
      ))
  }

  for (v in need_cols) {
    if (v %in% names(df)) next
    if (v %in% interaction_output_cols) next
    if (v %in% names(cov_vals)) {
      val     <- cov_vals[[v]]
      df[[v]] <- if (length(val) == n_obs) val else rep(val[[1L]], n_obs)
    } else {
      df[[v]] <- 0
      warning(sprintf(
        "build_cov_annual: column '%s' not in cov_vals-filled with 0.", v
      ))
    }
  }

  for (nm in names(interactions)) {
    cols    <- interactions[[nm]]
    df[[nm]] <- df[[cols[1L]]] * df[[cols[2L]]]
  }

  df
}


#' Extract names of non-seasonal covariates from a fitted model
#'
#' Returns the covariate names that appear in the model with a non-negligible
#' coefficient — i.e. those that actually affect the parameters and need to
#' be supplied during marginalization.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param free_cols Character vector of column names to always exclude from
#'   the result (intercept and seasonality terms are always excluded).
#'   Default \code{NULL}.
#' @param tol Numeric. Coefficients smaller than this in absolute value are
#'   treated as zero (LASSO shrinkage). Default \code{1e-4}.
#'
#' @return Character vector of active covariate names, or \code{NULL} if
#'   the model has no active non-seasonal covariates.
#'
#' @export
active_covariates <- function(fit, free_cols = NULL, tol = 1e-4) {

  if (!inherits(fit, "nhpp_fit"))
    stop("active_covariates: `fit` must be an nhpp_fit object.")

  always_exclude <- c("(Intercept)", "cos1", "sen1", "cos2", "sen2", free_cols)

  dm  <- fit$dm
  par <- fit$par

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  beta_mu  <- par[seq_len(p_mu)]
  beta_sig <- par[p_mu + seq_len(p_sig)]
  beta_xi  <- par[p_mu + p_sig + seq_len(p_xi)]

  names(beta_mu)  <- colnames(dm$X_mu)
  names(beta_sig) <- colnames(dm$X_sigma)
  names(beta_xi)  <- colnames(dm$X_xi)

  # Collect all candidate column names across all three blocks
  cn_mu  <- colnames(dm$X_mu) [!colnames(dm$X_mu)  %in% always_exclude]
  cn_sig <- colnames(dm$X_sigma)[!colnames(dm$X_sigma) %in% always_exclude]
  cn_xi  <- colnames(dm$X_xi) [!colnames(dm$X_xi)  %in% always_exclude]

  # Keep only those with a non-negligible coefficient
  active <- character(0L)
  for (cn in cn_mu)  if (abs(beta_mu[cn])  > tol) active <- c(active, cn)
  for (cn in cn_sig) if (cn %in% names(beta_sig) &&
                         abs(beta_sig[cn]) > tol) active <- c(active, cn)
  for (cn in cn_xi)  if (cn %in% names(beta_xi)  &&
                         abs(beta_xi[cn])  > tol) active <- c(active, cn)

  active <- unique(active)
  if (length(active) == 0L) NULL else active
}
