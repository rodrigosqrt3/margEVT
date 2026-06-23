# =============================================================================
# backtest.R
# Walk-forward backtesting for nhpp_fit objects.
#
# The validation windows are fully data-driven: the function finds the year
# range in `data`, starts validating after `min_train_years`, and steps
# forward by `window_years`. No hardcoded dates anywhere.
#
# For each window:
#   1. Fit the model on training years at the original lambda/alpha
#   2. Build empirical block bootstrap covariate sample from training years
#   3. Compute return levels on training data
#   4. Validate against observed annual maxima in validation years
#   5. Compute PIT (probability integral transform) for calibration
# =============================================================================

#' Walk-forward backtesting for NHPP extreme value models
#'
#' Evaluates predictive performance of a fitted \code{nhpp_fit} model via
#' walk-forward cross-validation. For each window, the model is refit on
#' training data at the original penalty and lambda, return levels are
#' computed, and exceedance rates are compared to nominal rates via
#' binomial tests. Probabilistic calibration is assessed via a
#' Kolmogorov-Smirnov test on the PIT values.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data The data frame used to fit \code{fit}. Must contain columns
#'   \code{y}, \code{year} (or \code{year_col}), and all active covariates.
#' @param varname Character. Name of the column in \code{data} containing
#'   the raw (non-declustered) series, used to compute observed annual maxima
#'   for validation.
#' @param TRs Numeric vector of return periods to evaluate. Default
#'   \code{c(2, 5, 10)}.
#' @param n_obs Integer. Observations per year. Default \code{365L}.
#' @param window_years Integer. Width of each validation window in years.
#'   Default \code{5L}.
#' @param min_train_years Integer. Minimum number of training years before
#'   the first validation window. Default \code{10L}.
#' @param n_boot Integer. Bootstrap years for covariate marginalization
#'   within each window. Default \code{200L}.
#' @param interactions Named list of interactions passed to
#'   \code{\link{build_cov_annual}}. Default \code{list()}.
#' @param year_col Character. Name of the year column. Default \code{"year"}.
#' @param min_obs_year Integer. Minimum observations in a validation year
#'   for it to be included. Default \code{as.integer(n_obs * 0.8)}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list with three elements:
#'   \item{results}{Data frame with one row per validation year, containing
#'     the observed annual maximum, PIT value, and exceedance indicators for
#'     each return period.}
#'   \item{binom_tests}{Data frame with binomial test results for each TR.}
#'   \item{ks_test}{Data frame with KS test result for PIT calibration,
#'     or \code{NULL} if fewer than 5 validation years had events.}
#'
#' @export
backtest <- function(fit, data, varname,
                     TRs             = c(2, 5, 10),
                     n_obs           = 365L,
                     window_years    = 5L,
                     min_train_years = 10L,
                     n_boot          = 200L,
                     interactions    = list(),
                     year_col        = "year",
                     min_obs_year    = NULL,
                     verbose         = TRUE) {

  if (!inherits(fit, "nhpp_fit"))
    stop("backtest: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("backtest: `data` must be a data frame.")
  if (!varname %in% names(data))
    stop(sprintf("backtest: column '%s' not found in data.", varname))
  if (!year_col %in% names(data))
    stop(sprintf("backtest: column '%s' not found in data.", year_col))
  if (!"y" %in% names(data))
    stop("backtest: `data` must contain column `y`.")

  n_obs        <- as.integer(n_obs)
  window_years <- as.integer(window_years)
  min_obs_year <- if (is.null(min_obs_year)) as.integer(n_obs * 0.8) else
    as.integer(min_obs_year)

  all_years  <- sort(unique(data[[year_col]]))
  min_year   <- min(all_years)
  max_year   <- max(all_years)
  thr        <- fit$threshold
  ac         <- active_covariates(fit)
  ac_in_data <- if (!is.null(ac)) ac[ac %in% names(data)] else character(0L)

  first_val <- min_year + min_train_years
  if (first_val > max_year)
    stop("backtest: not enough years for even one validation window given `min_train_years`.")

  val_starts <- seq(first_val, max_year, by = window_years)
  windows    <- lapply(val_starts, function(vs) {
    list(
      train_up_to = vs - 1L,
      val_start   = vs,
      val_end     = min(vs + window_years - 1L, max_year)
    )
  })

  results_list <- list()

  for (j in seq_along(windows)) {
    w         <- windows[[j]]
    train_yrs <- all_years[all_years <= w$train_up_to]
    val_yrs   <- all_years[all_years >= w$val_start &
                             all_years <= w$val_end]

    if (length(train_yrs) < 2L || length(val_yrs) == 0L) next

    if (verbose)
      message(sprintf("[Window %d] Train: %d-%d | Val: %d-%d",
                      j, min(train_yrs), max(train_yrs),
                      w$val_start, w$val_end))

    df_train <- data[data[[year_col]] %in% train_yrs, , drop = FALSE]

    fit_w <- .refit_boot(df_train, fit)
    if (is.null(fit_w) || !fit_w$converged) {
      if (verbose) message("  -> Did not converge, skipping window.")
      next
    }

    cov_w <- .build_cov_bootstrap(fit_w, data, train_yrs, ac_in_data,
                                  n_boot, n_obs, interactions, year_col)
    if (length(cov_w) == 0L) {
      if (verbose) message("  -> No valid covariate blocks, skipping window.")
      next
    }

    rl_w <- .compute_rl_from_cov_list(fit_w, cov_w, TRs, thr,
                                      max(data[[varname]], na.rm = TRUE) * 3,
                                      n_obs)

    if (verbose)
      message(sprintf("  RLs: %s",
                      paste(paste0("T", TRs, "=", round(rl_w, 1)),
                            collapse = " | ")))

    for (yr in val_yrs) {
      df_yr  <- data[data[[year_col]] == yr, , drop = FALSE]
      y_yr   <- df_yr[[varname]]
      exc_yr <- y_yr[!is.na(y_yr) & y_yr > thr]
      M_yr   <- if (length(exc_yr) > 0L) max(exc_yr) else NA_real_

      if (is.na(M_yr) || nrow(df_yr) < min_obs_year) {
        F_yr <- NA_real_
      } else {
        # PIT: P(annual max <= M_yr) under the training model
        F_yr <- mean(vapply(cov_w, function(cm)
          .f_annual_one(M_yr, fit_w, cm), numeric(1L)), na.rm = TRUE)
      }

      superou <- stats::setNames(
        vapply(rl_w, function(rl)
          if (!is.na(M_yr) && !is.na(rl)) M_yr > rl else FALSE,
          logical(1L)),
        paste0("exc_T", TRs)
      )

      row <- data.frame(
        window    = j,
        year      = yr,
        M_obs     = round(M_yr, 3L),
        F_ann     = round(F_yr, 4L),
        has_event = length(exc_yr) > 0L
      )
      row <- cbind(row, t(superou))
      results_list[[length(results_list) + 1L]] <- row
    }
  }

  if (length(results_list) == 0L) {
    warning("backtest: no validation results produced.")
    return(invisible(NULL))
  }

  df_results <- do.call(rbind, results_list)

  binom_list <- lapply(TRs, function(TR) {
    col   <- paste0("exc_T", TR)
    if (!col %in% names(df_results)) return(NULL)
    n_exc <- sum(df_results[[col]], na.rm = TRUE)
    n_tot <- sum(!is.na(df_results[[col]]))
    p_obs <- n_exc / n_tot
    p_teo <- 1 / TR
    bt    <- stats::binom.test(n_exc, n_tot, p = p_teo)
    data.frame(
      TR           = TR,
      expected     = round(p_teo, 4L),
      observed     = round(p_obs, 4L),
      n_exceedances = n_exc,
      n_total      = n_tot,
      p_value      = round(bt$p.value, 4L),
      calibrated   = bt$p.value > 0.05
    )
  })
  df_binom <- do.call(rbind, Filter(Negate(is.null), binom_list))

  if (verbose) {
    message("\n Binomial calibration tests")
    for (i in seq_len(nrow(df_binom))) {
      message(sprintf("  T=%3d | expected=%.3f | observed=%.3f (%d/%d) | p=%.3f %s",
                      df_binom$TR[i], df_binom$expected[i],
                      df_binom$observed[i], df_binom$n_exceedances[i],
                      df_binom$n_total[i], df_binom$p_value[i],
                      if (df_binom$calibrated[i]) "v" else "!"))
    }
  }

  df_ev  <- df_results[df_results$has_event & !is.na(df_results$F_ann), ]
  ks_res <- NULL
  if (nrow(df_ev) >= 5L) {
    ks    <- stats::ks.test(jitter(df_ev$F_ann, factor = 1e-6), "punif")
    ks_res <- data.frame(
      n           = nrow(df_ev),
      D_statistic = round(as.numeric(ks$statistic), 4L),
      p_value     = round(ks$p.value, 4L),
      calibrated  = ks$p.value > 0.05
    )
    if (verbose)
      message(sprintf(
        "\n PIT calibration (KS test) \n  n=%d | D=%.4f | p=%.4f | %s",
        ks_res$n, ks_res$D_statistic, ks_res$p_value,
        if (ks_res$calibrated) "Well calibrated" else "Miscalibration detected"
      ))
  }

  invisible(list(
    results     = df_results,
    binom_tests = df_binom,
    ks_test     = ks_res
  ))
}


# =============================================================================
# Internal helpers
# =============================================================================

# .build_cov_bootstrap()
# Sample n_boot years of covariates from training years.
.build_cov_bootstrap <- function(fit_w, data, train_yrs, ac_in_data,
                                 n_boot, n_obs, interactions, year_col) {

  if (length(ac_in_data) == 0L) {
    return(lapply(seq_len(n_boot), function(i)
      build_cov_annual(fit_w, list(), n_obs = n_obs,
                       interactions = interactions)))
  }

  col_need  <- unique(c(year_col, ac_in_data))
  col_need  <- col_need[col_need %in% names(data)]
  df_cov    <- data[data[[year_col]] %in% train_yrs,
                    col_need, drop = FALSE]

  # Impute NAs with column medians
  for (v in ac_in_data) {
    if (v %in% names(df_cov)) {
      nas <- is.na(df_cov[[v]])
      if (any(nas))
        df_cov[[v]][nas] <- stats::median(df_cov[[v]], na.rm = TRUE)
    }
  }

  anos_disp <- unique(df_cov[[year_col]][!is.na(df_cov[[year_col]])])
  if (length(anos_disp) == 0L) return(list())

  cov_list <- lapply(seq_len(n_boot), function(i) {
    yr    <- sample(anos_disp, 1L)
    df_yr <- df_cov[df_cov[[year_col]] == yr, ac_in_data, drop = FALSE]
    df_yr <- df_yr[stats::complete.cases(df_yr), , drop = FALSE]
    df_yr <- utils::head(df_yr, n_obs)
    if (nrow(df_yr) < as.integer(n_obs * 0.9)) return(NULL)
    if (nrow(df_yr) < n_obs) {
      falta <- n_obs - nrow(df_yr)
      df_yr <- rbind(df_yr,
                     df_yr[rep(nrow(df_yr), falta), , drop = FALSE])
    }
    build_cov_annual(fit_w, as.list(df_yr),
                     n_obs = n_obs, interactions = interactions)
  })
  Filter(Negate(is.null), cov_list)
}

# .f_annual_one()
# P(annual max <= z) for one covariate data frame.
.f_annual_one <- function(z, fit_w, cov_df) {
  params  <- predict_params(fit_w, as.data.frame(cov_df))
  mu_t    <- params$mu
  sigma_t <- params$sigma
  xi_t    <- params$xi
  xi_tol  <- 1e-6
  z_u     <- 1 + xi_t * (z - mu_t) / sigma_t
  lam_vec <- ifelse(
    abs(xi_t) < xi_tol,
    exp(pmax(-500, -(z - mu_t) / sigma_t)),
    ifelse(z_u <= 0, 0,
           exp((-1 / xi_t) * log(pmax(z_u, 1e-300))))
  )
  exp(-mean(lam_vec, na.rm = TRUE))
}

# .compute_rl_from_cov_list()
# Compute return levels by averaging over a list of covariate data frames.
.compute_rl_from_cov_list <- function(fit_w, cov_list, TRs, z_lo, z_hi,
                                      n_obs) {
  big_cov  <- as.data.frame(do.call(rbind, cov_list))
  params   <- predict_params(fit_w, big_cov)
  f_annual <- function(z)
    .annual_exceedance_prob(z, params$mu, params$sigma, params$xi, n_obs)
  sapply(TRs, function(TR)
    .find_return_level(TR, f_annual, z_lo, z_hi, tol = 1e-3))
}
