# =============================================================================
# marginalize.R
# Unconditional return level computation via marginalization over covariates.
#
# Three approaches:
#   A. Conditional   : fix covariates at user-specified scenarios
#   B. Parametric    : average over a user-supplied simulated covariate sample
#   C. Empirical     : resample observed complete covariate years
#
# For B and C the Monte Carlo standard error of each return level is
# reported, se = s_T / (sqrt(n) * g_T), where s_T is the standard deviation
# of the pathwise annual probabilities at the root and g_T the mean pathwise
# derivative with respect to the level.
# =============================================================================

# -----------------------------------------------------------------------------
# .tail_measure_at_level()
# Evaluates [1 + xi * (z - mu) / sigma]^(-1/xi), with the correct endpoint
# convention. Below the lower endpoint for xi > 0 the measure is infinite;
# above the upper endpoint for xi < 0 it is zero.
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

# Derivative of the tail measure with respect to the level z (non-positive).
# Zero above a finite upper endpoint; below a lower endpoint the annual
# probability is zero, so the contribution is irrelevant and set to zero.
.tail_measure_derivative <- function(z, mu_t, sigma_t, xi_t, xi_tol = 1e-6) {
  inner <- 1 + xi_t * (z - mu_t) / sigma_t
  out   <- numeric(length(inner))
  gumbel  <- abs(xi_t) < xi_tol
  regular <- !gumbel & inner > 0
  out[gumbel] <- -exp(pmax(-500, -(z - mu_t[gumbel]) / sigma_t[gumbel])) /
    sigma_t[gumbel]
  out[regular] <- -exp((-1 / xi_t[regular] - 1) * log(inner[regular])) /
    sigma_t[regular]
  out
}

.validate_param_stack <- function(mu_t, sigma_t, xi_t, n_obs) {
  n_obs <- as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L ||
      length(mu_t) != length(sigma_t) || length(mu_t) != length(xi_t) ||
      length(mu_t) == 0L || length(mu_t) %% n_obs != 0L)
    stop(".annual_exceedance_prob: parameter vectors must have equal, positive lengths divisible by `n_obs`.")
  if (any(!is.finite(mu_t)) || any(!is.finite(sigma_t)) ||
      any(sigma_t <= 0) || any(!is.finite(xi_t)))
    stop(".annual_exceedance_prob: model parameters must be finite and scales strictly positive.")
  n_obs
}

# Pathwise annual non-exceedance probabilities G(z | v_r) for a stack of
# n_mc trajectories, each with n_obs rows.
.pathwise_prob <- function(z, mu_t, sigma_t, xi_t, n_obs) {
  n_obs   <- .validate_param_stack(mu_t, sigma_t, xi_t, n_obs)
  lam_vec <- .tail_measure_at_level(z, mu_t, sigma_t, xi_t)
  if (any(is.na(lam_vec)) || any(lam_vec < 0))
    stop(".annual_exceedance_prob: invalid tail measure encountered.")
  lam_mat <- matrix(lam_vec, nrow = n_obs)
  exp(-colSums(lam_mat) / n_obs)
}

# Marginal annual non-exceedance probability G_Pi(z).
.annual_exceedance_prob <- function(z, mu_t, sigma_t, xi_t, n_obs) {
  mean(.pathwise_prob(z, mu_t, sigma_t, xi_t, n_obs))
}

# Conditional Monte Carlo standard error of the root z (Prop. in Villa and
# Ziegelmann, stochastic marginalization): s_T / (sqrt(n) g_T).
.mc_root_se <- function(z, mu_t, sigma_t, xi_t, n_obs) {
  if (!is.finite(z)) return(NA_real_)
  G  <- .pathwise_prob(z, mu_t, sigma_t, xi_t, n_obs)
  n  <- length(G)
  if (n < 2L) return(NA_real_)
  dl <- matrix(.tail_measure_derivative(z, mu_t, sigma_t, xi_t), nrow = n_obs)
  dG <- ifelse(G > 0, -G * colSums(dl) / n_obs, 0)
  g  <- mean(dG)
  if (!is.finite(g) || g <= 0) return(NA_real_)
  stats::sd(G) / (sqrt(n) * g)
}

# -----------------------------------------------------------------------------
# .find_return_level()
# Root-find z such that P(annual max <= z) = 1 - 1/TR.
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

# Return levels and Monte Carlo standard errors for one parameter stack.
.rl_from_params <- function(params, TRs, z_lo, z_hi, tol, n_obs,
                            with_se = TRUE) {
  f <- function(z)
    .annual_exceedance_prob(z, params$mu, params$sigma, params$xi, n_obs)
  rl <- vapply(TRs, function(TR) .find_return_level(TR, f, z_lo, z_hi, tol),
               numeric(1L))
  se <- if (with_se) vapply(rl, .mc_root_se, numeric(1L),
                            mu_t = params$mu, sigma_t = params$sigma,
                            xi_t = params$xi, n_obs = n_obs) else
    rep(NA_real_, length(TRs))
  list(rl = rl, se = se)
}

.validate_fixed <- function(fixed, fun = "marginalize") {
  if (!is.list(fixed) ||
      (length(fixed) > 0L &&
       (is.null(names(fixed)) || any(names(fixed) == "") ||
        anyDuplicated(names(fixed)) ||
        !all(vapply(fixed, function(v) is.numeric(v) && length(v) == 1L &&
                      is.finite(v), logical(1L))))))
    stop(sprintf("%s: `fixed` must be a named list of finite numeric scalars.",
                 fun))
  fixed
}


#' Compute unconditional return levels by marginalizing over covariates
#'
#' Given a fitted \code{nhpp_fit} model, computes return levels by
#' integrating out covariate uncertainty under one or more approaches:
#'
#' \describe{
#'   \item{A - conditional}{Covariates fixed at user-specified scenarios
#'     (e.g. quantiles of observed distributions). Fast; useful for
#'     sensitivity analysis.}
#'   \item{B - parametric}{Covariates drawn from a user-supplied Monte Carlo
#'     sample of complete annual trajectories (e.g. from
#'     \code{\link{simulate_covariates}}).}
#'   \item{C - empirical}{Annual-block resampling of observed covariate
#'     years from \code{data}, preserving within-year temporal structure.}
#' }
#'
#' For each approach and return period \eqn{T}, the return level solves
#' \eqn{G(z) = 1 - 1/T}, where \eqn{G} is the annual non-exceedance
#' probability averaged over trajectories. For approaches B and C the
#' conditional Monte Carlo standard error of the root,
#' \eqn{\hat s_T / (\sqrt{n}\,\hat g_T)}, is reported in \code{se_mc}, with
#' \eqn{\hat s_T} the standard deviation of the pathwise probabilities at the
#' root and \eqn{\hat g_T} the mean pathwise derivative with respect to the
#' level. It measures only finite-ensemble integration error, conditional on
#' the fitted model and trajectory law.
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
#'   \code{\link{build_cov_annual}}. If \code{NULL}, a single scenario named
#'   \code{"default"} with all active covariates set to zero is used.
#' @param mc_sample List of data frames for approach B. Each element is one
#'   simulated year of covariates (\code{n_obs} rows). Required if
#'   \code{"B" \%in\% approaches}.
#' @param interactions Named list of interactions passed to
#'   \code{\link{build_cov_annual}}. Default \code{list()}.
#' @param fixed Named list of scalar covariate values imposed on every
#'   scenario and trajectory, overriding values from \code{scenarios},
#'   \code{mc_sample}, and \code{data} (e.g. \code{list(trend = 1)} to
#'   evaluate return levels at the end of the record). Default \code{list()}.
#' @param year_col Name of the year column in \code{data}. Default
#'   \code{"year"}.
#' @param n_boot Integer. Number of resampled years for approach C.
#'   Default \code{500L}.
#' @param seed Integer seed for approach C (the caller's random-number stream
#'   is preserved). Default \code{2024L}.
#' @param tol_root Numeric. Root-finding tolerance. Default \code{1e-4}.
#' @param z_hi Numeric. Upper search bound for return levels. If \code{NULL},
#'   set to \code{3 * max(data$y, na.rm = TRUE)}.
#'
#' @return A data frame of class \code{nhpp_rl} with columns
#'   \code{approach}, \code{scenario}, \code{TR}, \code{RL}, and
#'   \code{se_mc} (\code{NA} for approach A). See \code{\link{plot.nhpp_rl}}
#'   and \code{\link{rl_table}}.
#'
#' @examples
#' set.seed(1)
#' n_years <- 30; n_obs <- 50
#' df <- data.frame(x = rnorm(n_years * n_obs),
#'                  year = rep(seq_len(n_years), each = n_obs))
#' df$y <- 10 + 2 * df$x + rexp(nrow(df), 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = n_obs, verbose = FALSE)
#' rl <- marginalize(fit, df, TRs = c(5, 10, 50),
#'                   scenarios = list(low = list(x = -1), high = list(x = 1)),
#'                   n_boot = 200)
#' rl
#' rl_table(rl)
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
                        fixed        = list(),
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
  fixed <- .validate_fixed(fixed)

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
  interaction_outputs <- names(interactions)
  relevant_interactions <- if (is.null(ac)) list() else
    interactions[interaction_outputs %in% ac]
  interaction_inputs <- unique(unlist(relevant_interactions,
                                      use.names = FALSE))
  direct_active <- if (is.null(ac)) character(0L) else
    setdiff(ac, interaction_outputs)
  required_covariates <- setdiff(unique(c(direct_active, interaction_inputs)),
                                 names(fixed))

  # Design columns that are inactive (zero operational coefficient) do not
  # affect the predicted parameters; they are filled with zero silently.
  design_cols <- unique(c(colnames(fit$dm$X_mu), colnames(fit$dm$X_sigma),
                          colnames(fit$dm$X_xi)))
  inactive <- setdiff(design_cols[design_cols != "(Intercept)" &
                                    !.is_harmonic(design_cols)],
                      c(required_covariates, interaction_outputs, names(fixed)))
  zero_fill <- stats::setNames(as.list(rep(0, length(inactive))), inactive)
  with_fixed <- function(vals)
    utils::modifyList(utils::modifyList(zero_fill, as.list(vals)), fixed)
  build <- function(vals)
    build_cov_annual(fit, with_fixed(vals), n_obs = n_obs, period = period,
                     interactions = interactions)
  out_rows <- function(app, scenario, res)
    data.frame(approach = app, scenario = scenario, TR = TRs,
               RL = round(res$rl, 4L), se_mc = res$se,
               stringsAsFactors = FALSE)

  results <- list()

  # A. Conditional scenarios
  if ("A" %in% approaches) {
    if (is.null(scenarios))
      scenarios <- list(default = stats::setNames(
        as.list(rep(0, length(required_covariates))), required_covariates))

    results[["A"]] <- do.call(rbind, lapply(names(scenarios), function(sc) {
      params <- predict_params(fit, build(scenarios[[sc]]))
      out_rows("A", sc, .rl_from_params(params, TRs, z_lo, z_hi, tol_root,
                                        n_obs, with_se = FALSE))
    }))
  }

  # B. Parametric Monte Carlo
  if ("B" %in% approaches) {
    missing_by_year <- vapply(mc_sample, function(yr_df)
      length(setdiff(required_covariates, names(yr_df))) > 0L, logical(1L))
    if (any(missing_by_year)) {
      first_bad <- which(missing_by_year)[1L]
      missing_names <- setdiff(required_covariates,
                               names(mc_sample[[first_bad]]))
      stop(sprintf(
        "marginalize: mc_sample[[%d]] is missing required covariates: %s",
        first_bad, paste(missing_names, collapse = ", ")
      ))
    }
    big_cov_B <- do.call(rbind, lapply(mc_sample, build))
    results[["B"]] <- out_rows(
      "B", "parametric_mc",
      .rl_from_params(predict_params(fit, big_cov_B), TRs, z_lo, z_hi,
                      tol_root, n_obs))
  }

  # C. Empirical annual-block resampling
  if ("C" %in% approaches) {
    if (!year_col %in% names(data))
      stop(sprintf("marginalize: column '%s' not found in data.", year_col))

    missing_required <- setdiff(required_covariates, names(data))
    if (length(missing_required) > 0L)
      stop("marginalize: approach C requires columns in `data`: ",
           paste(missing_required, collapse = ", "))

    if (length(required_covariates) == 0L) {
      # No time-varying covariates to resample: every year is identical.
      cov_C <- rep(list(build(list())), n_boot)
    } else {
      col_need   <- unique(c(year_col, required_covariates))
      df_valid   <- data[stats::complete.cases(data[, col_need, drop = FALSE]),
                         col_need, drop = FALSE]
      years_hist <- unique(df_valid[[year_col]])

      cov_C <- .with_preserved_seed(seed, {
        lapply(seq_len(n_boot), function(i) {
          yr    <- sample(years_hist, 1L)
          df_yr <- df_valid[df_valid[[year_col]] == yr,
                            required_covariates, drop = FALSE]
          df_yr <- utils::head(df_yr, n_obs)
          if (nrow(df_yr) < n_obs) return(NULL)
          build(df_yr)
        })
      })
      cov_C <- Filter(Negate(is.null), cov_C)

      if (length(cov_C) == 0L)
        stop("marginalize: approach C produced no valid bootstrap years.")
    }

    results[["C"]] <- out_rows(
      "C", "empirical_bootstrap",
      .rl_from_params(predict_params(fit, do.call(rbind, cov_C)), TRs, z_lo,
                      z_hi, tol_root, n_obs))
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  class(out) <- c("nhpp_rl", "data.frame")
  out
}


#' Plot return-level curves from marginalize()
#'
#' Draws return level against return period (logarithmic axis), one curve
#' per approach and scenario. For approaches with a Monte Carlo standard
#' error, dashed lines show \eqn{\pm 2} standard errors.
#'
#' @param x An \code{nhpp_rl} object returned by \code{\link{marginalize}}.
#' @param ... Additional graphical parameters passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{NULL}, called for side effects.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500), year = rep(1:30, each = 50))
#' df$y <- 10 + 2 * df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' plot(marginalize(fit, df, TRs = c(2, 5, 10, 50, 100), n_boot = 200))
#'
#' @export
plot.nhpp_rl <- function(x, ...) {
  x <- as.data.frame(x)
  keys <- unique(x[, c("approach", "scenario")])
  cols <- grDevices::hcl.colors(max(2L, nrow(keys)), "Dark 3")
  lo <- x$RL - 2 * ifelse(is.na(x$se_mc), 0, x$se_mc)
  hi <- x$RL + 2 * ifelse(is.na(x$se_mc), 0, x$se_mc)
  graphics::plot(range(x$TR), range(c(lo, hi), na.rm = TRUE), type = "n",
                 log = "x", xlab = "Return period (years)",
                 ylab = "Return level", main = "Return-level curves", ...)
  graphics::grid()
  labels <- character(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sub <- x[x$approach == keys$approach[i] & x$scenario == keys$scenario[i], ]
    sub <- sub[order(sub$TR), ]
    graphics::lines(sub$TR, sub$RL, type = "b", pch = 19, col = cols[i])
    if (any(!is.na(sub$se_mc))) {
      graphics::lines(sub$TR, sub$RL - 2 * sub$se_mc, lty = 2, col = cols[i])
      graphics::lines(sub$TR, sub$RL + 2 * sub$se_mc, lty = 2, col = cols[i])
    }
    labels[i] <- paste(keys$approach[i], keys$scenario[i], sep = ": ")
  }
  graphics::legend("topleft", legend = labels, col = cols[seq_len(nrow(keys))],
                   lty = 1, pch = 19, bty = "n")
  invisible(NULL)
}
