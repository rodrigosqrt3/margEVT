# =============================================================================
# marginalize.R
# Unconditional return level computation via marginalization over covariates.
#
# Three approaches:
#   A. Conditional   : fix covariates at user-specified scenarios
#   B. Parametric    : average over a user-supplied simulated covariate sample
#   C. Empirical     : block bootstrap over observed covariate years
#
# No climate-specific assumptions anywhere. The user supplies covariate
# scenarios, simulated trajectories, or historical data — this function
# just marginalizes over whatever it receives.
# =============================================================================

# -----------------------------------------------------------------------------
# .tail_measure_at_level()
# Evaluates [1 + xi * (z - mu) / sigma]^(-1/xi), with the correct endpoint
# convention. Below the lower endpoint for xi > 0 the measure is infinite;
# above the upper endpoint for xi < 0 it is zero.
# Internal — not exported.
# -----------------------------------------------------------------------------
.tail_measure_at_level <- function(z, mu_t, sigma_t, xi_t, xi_tol = 1e-6) {
  inner <- 1 + xi_t * (z - mu_t) / sigma_t
  out   <- rep(NA_real_, length(inner))

  gumbel <- is.finite(xi_t) & abs(xi_t) < xi_tol
  out[gumbel] <- exp(pmax(-500, -(z - mu_t[gumbel]) / sigma_t[gumbel]))

  regular <- is.finite(inner) & is.finite(xi_t) & !gumbel & inner > 0
  out[regular] <- exp((-1 / xi_t[regular]) * log(inner[regular]))

  below_lower <- is.finite(inner) & !gumbel & inner <= 0 & xi_t > 0
  above_upper <- is.finite(inner) & !gumbel & inner <= 0 & xi_t < 0
  out[below_lower] <- Inf
  out[above_upper] <- 0

  out
}

# -----------------------------------------------------------------------------
# .annual_exceedance_prob()
# Computes P(annual maximum <= z) for one stacked covariate matrix.
# Vectorized over n_mc Monte Carlo years stacked row-wise.
# Internal — not exported.
# -----------------------------------------------------------------------------
.annual_exceedance_prob <- function(z, mu_t, sigma_t, xi_t, n_obs) {
  n_obs <- as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L ||
      length(mu_t) != length(sigma_t) || length(mu_t) != length(xi_t) ||
      length(mu_t) == 0L || length(mu_t) %% n_obs != 0L)
    stop(".annual_exceedance_prob: parameter vectors must have equal, positive lengths divisible by `n_obs`.")
  lam_vec <- .tail_measure_at_level(z, mu_t, sigma_t, xi_t)

  n_mc    <- as.integer(length(mu_t) / n_obs)
  lam_mat <- matrix(lam_vec, nrow = n_obs, ncol = n_mc)

  mean(exp(-colSums(lam_mat, na.rm = TRUE) / n_obs), na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# .find_return_level()
# Root-find z such that P(annual max <= z) = 1 - 1/TR.
# Internal — not exported.
# -----------------------------------------------------------------------------
.find_return_level <- function(TR, f_annual, z_lo, z_hi, tol = 1e-4) {
  p_target <- 1 - 1 / TR
  objective <- function(z) f_annual(z) - p_target

  f_lo <- tryCatch(objective(z_lo), error = function(e) NA_real_)
  f_hi <- tryCatch(objective(z_hi), error = function(e) NA_real_)
  if (!is.finite(f_lo) || !is.finite(f_hi)) return(NA_real_)

  # The NHPP tail model is used only at and above its fitting threshold.
  # Do not let numerical bracketing extrapolate into the unidentified body.
  if (f_lo > 0) return(NA_real_)
  if (f_lo == 0) return(z_lo)

  # Expand only the upper endpoint when needed.
  span <- max(z_hi - z_lo, 1)
  n_expand <- 0L
  while (f_hi < 0 && n_expand < 20L) {
    n_expand <- n_expand + 1L
    z_hi <- z_hi + span * 2^(n_expand - 1L)
    f_hi <- tryCatch(objective(z_hi), error = function(e) NA_real_)
    if (!is.finite(f_hi)) return(NA_real_)
  }
  if (f_hi < 0) return(NA_real_)

  tryCatch(
    stats::uniroot(objective, lower = z_lo, upper = z_hi, tol = tol)$root,
    error = function(e) NA_real_
  )
}


#' Compute unconditional return levels by marginalizing over covariates
#'
#' Given a fitted \code{nhpp_fit} model, computes return levels by
#' integrating out covariate uncertainty under one or more approaches:
#'
#' \describe{
#'   \item{A — conditional}{Covariates fixed at user-specified scenarios
#'     (e.g. quantiles of observed distributions). Fast; useful for
#'     sensitivity analysis.}
#'   \item{B — parametric}{Covariates drawn from a user-supplied Monte Carlo
#'     sample (e.g. from a fitted VAR or any stochastic generator). The
#'     sample is a list of \code{n_obs}-row data frames, one per simulated
#'     year.}
#'   \item{C — empirical}{Block bootstrap over observed covariate years.
#'     Resamples whole years from \code{data} to preserve temporal structure.}
#' }
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data The data frame used to fit \code{fit}. Required for approach C
#'   and for computing \code{z_hi}. Must contain a column \code{year} (or
#'   the name given in \code{year_col}) and all active covariates.
#' @param TRs Numeric vector of return periods in years. Default
#'   \code{c(2, 5, 10, 20, 50, 100)}.
#' @param n_obs Integer observations per year. If \code{NULL}, uses the
#'   rounded \code{fit$obs_per_year} value.
#' @param period Numeric seasonal period in observations. If \code{NULL},
#'   inherits \code{fit$obs_per_year}.
#' @param approaches Character vector. Which approaches to run. Any subset of
#'   \code{c("A", "B", "C")}. Default \code{c("A", "C")}.
#' @param scenarios Named list of covariate scenarios for approach A. Each
#'   element is itself a named list of scalar covariate values, passed to
#'   \code{\link{build_cov_annual}}. If \code{NULL}, a single scenario with
#'   all covariates set to zero is used.
#' @param mc_sample List of data frames for approach B. Each element is one
#'   simulated year of covariates (\code{n_obs} rows). Required if
#'   \code{"B" \%in\% approaches}.
#' @param interactions Named list of interactions passed to
#'   \code{\link{build_cov_annual}}. Default \code{list()}.
#' @param year_col Name of the year column in \code{data}. Default
#'   \code{"year"}.
#' @param n_boot Integer. Number of bootstrap years for approach C.
#'   Default \code{500L}.
#' @param seed Integer. Random seed for approach C. Default \code{2024L}.
#' @param tol_root Numeric. Root-finding tolerance. Default \code{1e-4}.
#' @param z_hi Numeric. Upper search bound for return levels. If \code{NULL},
#'   set to \code{3 * max(data$y, na.rm = TRUE)}.
#'
#' @return A data frame with columns \code{approach}, \code{scenario},
#'   \code{TR}, \code{RL}.
#'
#' @export
marginalize <- function(fit, data,
                        TRs          = c(2, 5, 10, 20, 50, 100),
                        n_obs        = NULL,
                        period       = NULL,
                        approaches   = c("A", "C"),
                        scenarios    = NULL,
                        mc_sample    = NULL,
                        interactions = list(),
                        year_col     = "year",
                        n_boot       = 500L,
                        seed         = 2024L,
                        tol_root     = 1e-4,
                        z_hi         = NULL) {

  if (!inherits(fit, "nhpp_fit"))
    stop("marginalize: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("marginalize: `data` must be a data frame.")
  if (!"y" %in% names(data) || !is.numeric(data$y) ||
      all(!is.finite(data$y)))
    stop("marginalize: `data$y` must contain at least one finite numeric value.")
  approaches <- unique(approaches)
  if (!is.character(approaches) || length(approaches) < 1L ||
      any(!approaches %in% c("A", "B", "C")))
    stop("marginalize: `approaches` must contain only 'A', 'B', or 'C'.")
  if (!is.numeric(TRs) || length(TRs) < 1L ||
      any(!is.finite(TRs)) || any(TRs <= 1))
    stop("marginalize: `TRs` must contain finite return periods greater than 1.")
  if ("B" %in% approaches && is.null(mc_sample))
    stop("marginalize: approach B requires `mc_sample` (list of covariate data frames).")
  if ("B" %in% approaches &&
      (!is.list(mc_sample) || length(mc_sample) < 1L ||
       !all(vapply(mc_sample, is.data.frame, logical(1L)))))
    stop("marginalize: `mc_sample` must be a non-empty list of data frames.")

  inherited_obs <- if (is.null(fit$obs_per_year)) 365.25 else fit$obs_per_year
  n_obs  <- if (is.null(n_obs)) as.integer(round(inherited_obs)) else
    as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L)
    stop("marginalize: `n_obs` must be a positive integer.")
  n_boot <- as.integer(n_boot)
  if (!is.finite(n_boot) || n_boot < 1L)
    stop("marginalize: `n_boot` must be a positive integer.")
  z_lo   <- fit$threshold
  z_hi   <- if (is.null(z_hi)) 3 * max(data$y, na.rm = TRUE) else z_hi
  if (!is.numeric(z_hi) || length(z_hi) != 1L || !is.finite(z_hi) ||
      z_hi <= z_lo)
    stop("marginalize: `z_hi` must be a finite value greater than the threshold.")
  ac     <- active_covariates(fit)

  results <- list()

  # A. Conditional scenarios
  if ("A" %in% approaches) {
    if (is.null(scenarios))
      scenarios <- list("default" = list())

    res_A <- do.call(rbind, lapply(names(scenarios), function(sc_name) {
      cov_df   <- build_cov_annual(fit, scenarios[[sc_name]],
                                   n_obs = n_obs,
                                   period = period,
                                   interactions = interactions)
      params   <- predict_params(fit, cov_df)
      f_annual <- function(z)
        .annual_exceedance_prob(z, params$mu, params$sigma, params$xi, n_obs)
      rl <- sapply(TRs, function(TR)
        .find_return_level(TR, f_annual, z_lo, z_hi, tol_root))
      data.frame(approach = "A", scenario = sc_name,
                 TR = TRs, RL = round(rl, 4L))
    }))
    results[["A"]] <- res_A
  }

  # B. Parametric Monte Carlo
  if ("B" %in% approaches) {
    cov_B     <- lapply(mc_sample, function(yr_df)
      build_cov_annual(fit, as.list(yr_df), n_obs = n_obs,
                       period = period,
                       interactions = interactions))
    big_cov_B <- as.data.frame(do.call(rbind, cov_B))
    params_B  <- predict_params(fit, big_cov_B)
    f_B       <- function(z)
      .annual_exceedance_prob(z, params_B$mu, params_B$sigma,
                              params_B$xi, n_obs)
    rl_B <- sapply(TRs, function(TR)
      .find_return_level(TR, f_B, z_lo, z_hi, tol_root))
    results[["B"]] <- data.frame(
      approach = "B", scenario = "parametric_mc",
      TR = TRs, RL = round(rl_B, 4L)
    )
  }

  # C. Empirical block bootstrap
  if ("C" %in% approaches) {
    if (!year_col %in% names(data))
      stop(sprintf("marginalize: column '%s' not found in data.", year_col))

    ac_in_data <- if (!is.null(ac)) ac[ac %in% names(data)] else character(0L)

    if (length(ac_in_data) == 0L) {
      # Stationary model: no covariates to resample, just build one empty frame
      cov_C <- lapply(seq_len(n_boot), function(i)
        build_cov_annual(fit, list(), n_obs = n_obs,
                         period = period,
                         interactions = interactions))
    } else {
      col_need   <- unique(c(year_col, ac_in_data))
      df_valid   <- data[stats::complete.cases(data[, col_need, drop = FALSE]),
                         col_need, drop = FALSE]
      years_hist <- unique(df_valid[[year_col]])

      cov_C <- .with_preserved_seed(seed, {
        lapply(seq_len(n_boot), function(i) {
          yr    <- sample(years_hist, 1L)
          df_yr <- df_valid[df_valid[[year_col]] == yr,
                            ac_in_data, drop = FALSE]
          df_yr <- utils::head(df_yr, n_obs)
          if (nrow(df_yr) < n_obs) return(NULL)
          build_cov_annual(fit, as.list(df_yr),
                           n_obs = n_obs, period = period,
                           interactions = interactions)
        })
      })
      cov_C <- Filter(Negate(is.null), cov_C)

      if (length(cov_C) == 0L)
        stop("marginalize: approach C produced no valid bootstrap years.")
    }

    big_cov_C <- as.data.frame(do.call(rbind, cov_C))
    params_C  <- predict_params(fit, big_cov_C)
    f_C       <- function(z)
      .annual_exceedance_prob(z, params_C$mu, params_C$sigma,
                              params_C$xi, n_obs)
    rl_C <- sapply(TRs, function(TR)
      .find_return_level(TR, f_C, z_lo, z_hi, tol_root))
    results[["C"]] <- data.frame(
      approach = "C", scenario = "empirical_bootstrap",
      TR = TRs, RL = round(rl_C, 4L)
    )
  }

  do.call(rbind, results)
}
