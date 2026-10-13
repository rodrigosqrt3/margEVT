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
#   5. Compute a tail-conditional PIT diagnostic
# =============================================================================

#' Walk-forward backtesting for NHPP extreme value models
#'
#' Evaluates predictive performance of a fitted \code{nhpp_fit} model via
#' walk-forward cross-validation. For each window, the model is refit on
#' training data at the original penalty and lambda, return levels are
#' computed, and exceedance rates are compared to nominal rates via
#' binomial tests. A Kolmogorov-Smirnov comparison on tail-conditional PIT
#' values is returned as an exploratory diagnostic; it is not a formal test
#' of global predictive calibration. Validation years with fewer than
#' \code{min_obs_year} observed values are excluded from binomial denominators.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data The data frame used to fit \code{fit}. Must contain columns
#'   \code{y}, \code{year} (or \code{year_col}), and all active covariates.
#' @param varname Character. Name of the column in \code{data} containing
#'   the raw (non-declustered) series, used to compute observed annual maxima
#'   for validation.
#' @param TRs Numeric vector of return periods to evaluate. Default
#'   \code{c(2, 5, 10)}.
#' @param n_obs Integer observations per year. If \code{NULL}, uses the
#'   rounded \code{fit$obs_per_year} value.
#' @param period Numeric seasonal period in observations. If \code{NULL},
#'   inherits \code{fit$obs_per_year}.
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
#' @param seed Integer. Seed used for empirical annual-block sampling.
#'   Default \code{2024L}. Each validation window uses a distinct derived seed.
#' @param trend_col Optional name of a linear time-trend covariate. Within
#'   each window it is recomputed as \eqn{(t - t_0)/(t_1 - t_0)}, where
#'   \eqn{t_0} is the first year of \code{data} and \eqn{t_1} the last
#'   training year, so that the trend equals 1 at the end of training.
#'   Return levels are evaluated at trend 1; the PIT of a validation year
#'   uses its extrapolated trend value. Default \code{NULL} (no trend).
#' @param fixed Named list of scalar covariate values imposed on every
#'   covariate block (see \code{\link{marginalize}}). Default \code{list()}.
#' @param reselect_lambda Logical. If \code{TRUE}, the penalty is re-selected
#'   by BIC in every training window; otherwise the original \code{lambda} is
#'   kept. Adaptive penalty weights are always recomputed within each training
#'   window, so the weights do not use validation data. Default \code{FALSE}.
#' @param cores Positive integer number of worker processes over which the
#'   windows are distributed. Results do not depend on \code{cores}.
#'   Default \code{1L}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list with three elements:
#'   \item{results}{Data frame with one row per validation year, containing
#'     the observed annual maximum, PIT value, and exceedance indicators for
#'     each return period.}
#'   \item{binom_tests}{Data frame with nominal binomial comparisons for each
#'     TR. The \code{not_rejected} field is the recommended interpretation;
#'     \code{calibrated} is retained as a legacy alias.}
#'   \item{ks_test}{Data frame with the exploratory tail-conditional PIT
#'     KS comparison, including \code{not_rejected} and the legacy
#'     \code{calibrated} alias, or \code{NULL} if fewer than 5 validation
#'     years had events.}
#'
#' @seealso \code{\link{marginalize}}, \code{\link{fit_nhpp}}
#'
#' @examples
#' set.seed(1)
#' n_years <- 30; n_obs <- 50
#' df <- data.frame(x = rnorm(n_years * n_obs),
#'                  year = rep(1991:2020, each = n_obs))
#' df$y <- 5 + df$x + rexp(nrow(df), 1 / 2)
#' fit <- fit_nhpp(df, threshold = 10, loc_vars = "x", penalty = "none",
#'                 obs_per_year = n_obs, verbose = FALSE)
#' \donttest{
#' bt <- backtest(fit, df, varname = "y", TRs = c(2, 5), n_obs = n_obs,
#'                window_years = 5, min_train_years = 15, n_boot = 50,
#'                verbose = FALSE)
#' bt$binom_tests
#' }
#'
#' @export
backtest <- function(fit, data, varname,
                     TRs             = c(2, 5, 10),
                     n_obs           = NULL,
                     period          = NULL,
                     window_years    = 5L,
                     min_train_years = 10L,
                     n_boot          = 200L,
                     interactions    = list(),
                     year_col        = "year",
                     min_obs_year    = NULL,
                     seed            = 2024L,
                     trend_col       = NULL,
                     fixed           = list(),
                     reselect_lambda = FALSE,
                     cores           = 1L,
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
  if (!is.numeric(TRs) || length(TRs) < 1L ||
      any(!is.finite(TRs)) || any(TRs <= 1))
    stop("backtest: `TRs` must contain finite return periods greater than 1.")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
    stop("backtest: `seed` must be a single finite numeric value.")
  if (!is.logical(reselect_lambda) || length(reselect_lambda) != 1L ||
      is.na(reselect_lambda))
    stop("backtest: `reselect_lambda` must be TRUE or FALSE.")
  cores <- .validate_cores(cores, "backtest")
  .validate_fixed(fixed, "backtest")
  model_cols <- unique(c(colnames(fit$dm$X_mu), colnames(fit$dm$X_sigma),
                         colnames(fit$dm$X_xi)))
  if (!is.null(trend_col) &&
      (!is.character(trend_col) || length(trend_col) != 1L ||
       !trend_col %in% model_cols))
    stop("backtest: `trend_col` must name a covariate of the fitted model.")

  inherited_obs <- if (is.null(fit$obs_per_year)) 365.25 else fit$obs_per_year
  n_obs        <- if (is.null(n_obs)) as.integer(round(inherited_obs)) else
    as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L)
    stop("backtest: `n_obs` must be a positive integer.")
  window_years <- as.integer(window_years)
  min_train_years <- as.integer(min_train_years)
  n_boot <- as.integer(n_boot)
  if (!is.finite(window_years) || window_years < 1L)
    stop("backtest: `window_years` must be a positive integer.")
  if (!is.finite(min_train_years) || min_train_years < 1L)
    stop("backtest: `min_train_years` must be a positive integer.")
  if (!is.finite(n_boot) || n_boot < 1L)
    stop("backtest: `n_boot` must be a positive integer.")
  min_obs_year <- if (is.null(min_obs_year)) as.integer(n_obs * 0.8) else
    as.integer(min_obs_year)
  if (!is.finite(min_obs_year) || min_obs_year < 1L)
    stop("backtest: `min_obs_year` must be a positive integer.")

  all_years  <- sort(unique(data[[year_col]][!is.na(data[[year_col]])]))
  if (length(all_years) == 0L)
    stop("backtest: `year_col` contains no valid years.")
  min_year   <- min(all_years)
  max_year   <- max(all_years)
  thr        <- fit$threshold
  ac         <- active_covariates(fit)
  interaction_outputs <- names(interactions)
  relevant_interactions <- if (is.null(ac)) list() else
    interactions[interaction_outputs %in% ac]
  interaction_inputs <- unique(unlist(relevant_interactions,
                                      use.names = FALSE))
  direct_active <- if (is.null(ac)) character(0L) else
    setdiff(ac, interaction_outputs)
  ac_in_data <- setdiff(unique(c(direct_active, interaction_inputs)),
                        c(trend_col, names(fixed)))
  missing_covariates <- setdiff(ac_in_data, names(data))
  if (length(missing_covariates) > 0L)
    stop("backtest: required covariates not found in data: ",
         paste(missing_covariates, collapse = ", "))

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

  # Adaptive weights are recomputed inside every training window (no look-ahead);
  # user-supplied weights are kept.
  weights_source <- if (is.null(fit$weights_source)) "none" else
    fit$weights_source
  window_weights <- if (startsWith(weights_source, "adaptive")) "adaptive" else
    NULL
  window_lambda <- if (reselect_lambda && fit$penalty != "none") "bic" else
    fit$lambda
  z_hi <- max(data[[varname]], na.rm = TRUE) * 3
  covariate_cols <- model_cols[model_cols != "(Intercept)" &
                                 !.is_harmonic(model_cols)]

  # Value of the trend covariate in `years` for a window whose training
  # period ends in `train_end`; equals 1 at the end of training.
  trend_value <- function(years, train_end)
    (years - min_year) / (train_end - min_year)

  # Replaces the trend column of a covariate block and refreshes the
  # interactions that use it.
  set_trend <- function(cm, value) {
    cm[[trend_col]] <- value
    for (nm in names(interactions))
      if (trend_col %in% interactions[[nm]])
        cm[[nm]] <- cm[[interactions[[nm]][1L]]] * cm[[interactions[[nm]][2L]]]
    cm
  }

  # One walk-forward window: returns validation rows and progress messages.
  # Deterministic given `j`, so the result does not depend on `cores`.
  run_window <- function(j) {
    w         <- windows[[j]]
    train_yrs <- all_years[all_years <= w$train_up_to]
    val_yrs   <- all_years[all_years >= w$val_start &
                             all_years <= w$val_end]
    msgs <- character(0L)
    if (length(train_yrs) < 2L || length(val_yrs) == 0L)
      return(list(rows = list(), msgs = msgs))

    msgs <- sprintf("[Window %d] Train: %d-%d | Val: %d-%d",
                    j, min(train_yrs), max(train_yrs), w$val_start, w$val_end)
    train_end <- max(train_yrs)

    df_train <- data[data[[year_col]] %in% train_yrs, , drop = FALSE]
    fixed_w  <- fixed
    if (!is.null(trend_col)) {
      df_train[[trend_col]] <- trend_value(df_train[[year_col]], train_end)
      fixed_w[[trend_col]]  <- 1
    }

    fit_w <- .refit_boot(df_train, fit, lambda = window_lambda,
                         penalty_weights = window_weights)
    if (is.null(fit_w) || !fit_w$converged)
      return(list(rows = list(),
                  msgs = c(msgs, "  -> Did not converge, skipping window.")))

    # Covariates resampled for this window: those required by the original
    # fit plus any that are active in the window refit. Columns inactive in
    # the window refit do not affect the parameters and are set to zero.
    ac_w <- active_covariates(fit_w)
    inputs_w <- unique(c(
      setdiff(ac_w, interaction_outputs),
      unlist(interactions[interaction_outputs %in% ac_w], use.names = FALSE)
    ))
    cols_w <- setdiff(unique(c(ac_in_data, intersect(inputs_w, names(data)))),
                      c(trend_col, names(fixed)))
    inactive_w <- setdiff(covariate_cols,
                          c(cols_w, interaction_outputs, trend_col, names(fixed)))
    fixed_w <- utils::modifyList(
      stats::setNames(as.list(rep(0, length(inactive_w))), inactive_w),
      fixed_w)

    cov_w <- .build_cov_bootstrap(
      fit_w        = fit_w,
      data         = data,
      train_yrs    = train_yrs,
      ac_in_data   = cols_w,
      n_boot       = n_boot,
      n_obs        = n_obs,
      interactions = interactions,
      year_col     = year_col,
      period       = period,
      seed         = as.integer(seed) + j - 1L,
      fixed        = fixed_w
    )
    if (length(cov_w) == 0L)
      return(list(rows = list(),
                  msgs = c(msgs, "  -> No valid covariate blocks, skipping window.")))

    rl_w <- .compute_rl_from_cov_list(fit_w, cov_w, TRs, thr, z_hi, n_obs)
    msgs <- c(msgs, sprintf("  RLs: %s",
                            paste(paste0("T", TRs, "=", round(rl_w, 1)),
                                  collapse = " | ")))

    rows <- lapply(val_yrs, function(yr) {
      df_yr  <- data[data[[year_col]] == yr, , drop = FALSE]
      y_yr   <- df_yr[[varname]]
      exc_yr <- y_yr[!is.na(y_yr) & y_yr > thr]
      M_yr   <- if (length(exc_yr) > 0L) max(exc_yr) else NA_real_
      valid_year <- sum(!is.na(y_yr)) >= min_obs_year

      if (is.na(M_yr) || !valid_year) {
        F_yr <- NA_real_
      } else {
        # PIT: P(annual max <= M_yr) under the training model, with the
        # trend extrapolated to the validation year.
        cov_yr <- if (is.null(trend_col)) cov_w else
          lapply(cov_w, set_trend, value = trend_value(yr, train_end))
        F_vals <- vapply(cov_yr, function(cm)
          .f_annual_one(M_yr, fit_w, cm), numeric(1L))
        F_yr <- if (all(is.finite(F_vals))) mean(F_vals) else NA_real_
      }

      superou <- stats::setNames(
        vapply(rl_w, function(rl)
          if (!valid_year || is.na(rl)) NA else
            if (is.na(M_yr)) FALSE else M_yr > rl,
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
      cbind(row, t(superou))
    })
    list(rows = rows, msgs = msgs)
  }

  window_out <- .par_lapply(seq_along(windows), run_window, cores = cores)
  if (verbose)
    for (wo in window_out) for (m in wo$msgs) message(m)
  results_list <- unlist(lapply(window_out, `[[`, "rows"), recursive = FALSE)

  if (length(results_list) == 0L) {
    warning("backtest: no validation results produced.")
    return(invisible(NULL))
  }

  df_results <- do.call(rbind, results_list)

  binom_list <- lapply(TRs, function(TR) {
    col   <- paste0("exc_T", TR)
    n_exc <- sum(df_results[[col]], na.rm = TRUE)
    n_tot <- sum(!is.na(df_results[[col]]))
    p_teo <- 1 / TR
    if (n_tot == 0L) {
      return(data.frame(
        TR = TR, expected = round(p_teo, 4L), observed = NA_real_,
        n_exceedances = 0L, n_total = 0L, p_value = NA_real_,
        not_rejected = NA, calibrated = NA
      ))
    }
    p_obs <- n_exc / n_tot
    bt    <- stats::binom.test(n_exc, n_tot, p = p_teo)
    data.frame(
      TR           = TR,
      expected     = round(p_teo, 4L),
      observed     = round(p_obs, 4L),
      n_exceedances = n_exc,
      n_total      = n_tot,
      p_value      = round(bt$p.value, 4L),
      not_rejected = bt$p.value > 0.05,
      calibrated   = bt$p.value > 0.05
    )
  })
  df_binom <- do.call(rbind, Filter(Negate(is.null), binom_list))

  if (verbose) {
    message("\n Nominal binomial rate comparisons")
    for (i in seq_len(nrow(df_binom))) {
      message(sprintf("  T=%3d | expected=%.3f | observed=%.3f (%d/%d) | p=%.3f %s",
                      df_binom$TR[i], df_binom$expected[i],
                      df_binom$observed[i], df_binom$n_exceedances[i],
                      df_binom$n_total[i], df_binom$p_value[i],
                      if (is.na(df_binom$not_rejected[i])) "not estimable" else
                        if (df_binom$not_rejected[i]) "compatible" else "departure"))
    }
  }

  df_ev  <- df_results[df_results$has_event & !is.na(df_results$F_ann), ]
  ks_res <- NULL
  if (nrow(df_ev) >= 5L) {
    pit_jittered <- .with_preserved_seed(
      as.integer(seed) + 100000L,
      jitter(df_ev$F_ann, factor = 1e-6)
    )
    ks <- stats::ks.test(pit_jittered, "punif")
    ks_res <- data.frame(
      n           = nrow(df_ev),
      D_statistic = round(as.numeric(ks$statistic), 4L),
      p_value     = round(ks$p.value, 4L),
      not_rejected = ks$p.value > 0.05,
      calibrated   = ks$p.value > 0.05
    )
    if (verbose)
      message(sprintf(
        "\n Exploratory tail-PIT KS comparison \n  n=%d | D=%.4f | p=%.4f | %s",
        ks_res$n, ks_res$D_statistic, ks_res$p_value,
        if (ks_res$not_rejected) "no nominal departure detected" else "nominal departure detected"
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
                                 n_boot, n_obs, interactions, year_col,
                                 period = NULL, seed = 2024L,
                                 fixed = list()) {

  if (length(ac_in_data) == 0L) {
    return(lapply(seq_len(n_boot), function(i)
      build_cov_annual(fit_w, fixed, n_obs = n_obs,
                       period = period,
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

  cov_list <- .with_preserved_seed(seed, {
    lapply(seq_len(n_boot), function(i) {
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
      build_cov_annual(fit_w, utils::modifyList(as.list(df_yr), fixed),
                       n_obs = n_obs, period = period,
                       interactions = interactions)
    })
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
  .annual_exceedance_prob(z, mu_t, sigma_t, xi_t,
                          n_obs = length(mu_t))
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
