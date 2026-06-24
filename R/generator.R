# =============================================================================
# generator.R
# Generic VAR-based stochastic generator over active covariates, for use
# as approach B (`mc_sample`) in marginalize().
#
# Workflow:
#   gen <- fit_var_generator(fit, data)
#   mc  <- simulate_covariates(gen, n_mc = 3000L, seed = 2024L)
#   marginalize(fit, data, approaches = c("A","B","C"), mc_sample = mc)
# =============================================================================

#' Fit a generic VAR generator over a model's active covariates
#'
#' Deseasonalizes each active covariate column in \code{data} (regressing
#' on \code{cos1, sen1, cos2, sen2}), then fits a multivariate VAR to the
#' standardized residuals. No covariate names are hardcoded: whatever
#' \code{\link{active_covariates}} returns (or whatever you pass in
#' \code{vars}) is treated as a final, ready-to-model column.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data Data frame used to fit \code{fit} (or any data frame
#'   containing the same covariate columns, plus seasonal harmonics if
#'   you already have them).
#' @param vars Character vector of covariate names to model. Default
#'   \code{NULL}: uses \code{active_covariates(fit)}, intersected with
#'   \code{names(data)}.
#' @param period Numeric. Seasonal period in observations. Default
#'   \code{365.25}.
#' @param lag_max Integer. Max VAR lag considered by \code{VARselect}.
#'   Default \code{NULL}: \code{min(5, floor(n / (3*K)))}.
#'
#' @return An object of class \code{nhpp_var_generator}, or \code{NULL}
#'   (with a warning) if the model has no active covariates.
#'
#' @export
fit_var_generator <- function(fit, data, vars = NULL, period = 365.25,
                              lag_max = NULL) {

  if (!inherits(fit, "nhpp_fit"))
    stop("fit_var_generator: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("fit_var_generator: `data` must be a data frame.")

  if (!requireNamespace("vars", quietly = TRUE))
    stop("fit_var_generator: package 'vars' is required.")

  if (is.null(vars)) {
    ac   <- active_covariates(fit)
    vars <- if (is.null(ac)) character(0L) else ac[ac %in% names(data)]
  } else {
    missing_v <- setdiff(vars, names(data))
    if (length(missing_v) > 0L)
      stop("fit_var_generator: columns not in data: ",
           paste(missing_v, collapse = ", "))
  }

  if (length(vars) == 0L) {
    warning("fit_var_generator: no active covariates - returning NULL generator.")
    return(NULL)
  }

  n <- nrow(data)
  if (!all(c("cos1", "sen1", "cos2", "sen2") %in% names(data))) {
    t_grid     <- seq_len(n)
    data$cos1  <- cos(2 * pi * t_grid / period)
    data$sen1  <- sin(2 * pi * t_grid / period)
    data$cos2  <- cos(4 * pi * t_grid / period)
    data$sen2  <- sin(4 * pi * t_grid / period)
  }

  cols_needed <- unique(c(vars, "cos1", "sen1", "cos2", "sen2"))
  df_base     <- data[, cols_needed, drop = FALSE]
  df_base     <- df_base[stats::complete.cases(df_base), ]

  if (nrow(df_base) < 10L)
    stop("fit_var_generator: too few complete rows after removing NAs.")

  seasonal_models <- list()
  resid_mat        <- matrix(NA_real_, nrow = nrow(df_base), ncol = length(vars))
  colnames(resid_mat) <- vars

  for (v in vars) {
    frm    <- stats::as.formula(paste0("`", v, "` ~ cos1 + sen1 + cos2 + sen2"))
    f      <- stats::lm(frm, data = df_base)
    res    <- stats::residuals(f)
    mu_res <- mean(res)
    sd_res <- stats::sd(res)
    if (!is.finite(sd_res) || sd_res < 1e-8) sd_res <- 1
    seasonal_models[[v]] <- list(coefs = stats::coef(f), mu_res = mu_res, sd_res = sd_res)
    resid_mat[, v] <- (res - mu_res) / sd_res
  }

  K <- length(vars)
  resid_sel <- resid_mat
  if (K == 1L) {
    # vars::VAR needs K >= 2; pad with noise, drop it again at simulation time
    resid_sel <- cbind(resid_mat, .dummy_noise = stats::rnorm(nrow(resid_mat)))
  }

  if (is.null(lag_max))
    lag_max <- max(1L, min(5L, floor(nrow(resid_sel) / (3L * ncol(resid_sel)))))

  p_opt <- as.integer(
    vars::VARselect(resid_sel, lag.max = lag_max, type = "none")$selection["SC(n)"]
  )
  p_opt <- max(1L, p_opt)

  fit_var <- vars::VAR(resid_sel, p = p_opt, type = "none")

  structure(
    list(
      vars            = vars,
      seasonal_models = seasonal_models,
      fit_var         = fit_var,
      p_opt           = p_opt,
      period          = period,
      var_colnames    = colnames(resid_sel)
    ),
    class = "nhpp_var_generator"
  )
}


#' Simulate stationary covariate trajectories from a VAR generator
#'
#' @param generator Object from \code{\link{fit_var_generator}}. If
#'   \code{NULL}, returns \code{NULL} (stationary/no-covariate model).
#' @param n_mc Integer. Number of simulated years.
#' @param n_obs Integer. Observations per simulated year. Default \code{365L}.
#' @param burn_in Integer. Burn-in to discard transient VAR dynamics.
#'   Default \code{300L}.
#' @param seed Integer. Random seed. Default \code{NULL} (not set).
#'
#' @return A list of length \code{n_mc}, each element an \code{n_obs}-row
#'   data frame with one column per active covariate - directly usable as
#'   \code{mc_sample} in \code{\link{marginalize}} (approach B).
#'
#' @export
simulate_covariates <- function(generator, n_mc, n_obs = 365L,
                                burn_in = 300L, seed = NULL) {

  if (is.null(generator)) return(NULL)
  if (!inherits(generator, "nhpp_var_generator"))
    stop("simulate_covariates: `generator` must come from fit_var_generator().")

  if (!is.null(seed)) set.seed(seed)

  n_obs  <- as.integer(n_obs)
  n_days <- n_mc * n_obs
  K      <- generator$fit_var$K
  p      <- generator$fit_var$p
  A      <- vars::Bcoef(generator$fit_var)
  resids <- stats::residuals(generator$fit_var)
  Sigma  <- crossprod(resids) / nrow(resids)
  chol_S <- t(chol(Sigma + diag(1e-8, K)))

  sim_pad <- matrix(0, nrow = n_days + burn_in + p, ncol = K)
  colnames(sim_pad) <- generator$var_colnames

  for (t in (p + 1L):nrow(sim_pad)) {
    y_lag <- numeric(K * p)
    for (i in seq_len(p))
      y_lag[((i - 1L) * K + 1L):(i * K)] <- sim_pad[t - i, ]
    sim_pad[t, ] <- as.numeric(A %*% y_lag) + as.numeric(chol_S %*% stats::rnorm(K))
  }
  sim_pad <- sim_pad[(burn_in + p + 1L):nrow(sim_pad), , drop = FALSE]

  day_grid <- rep(seq_len(n_obs), times = n_mc)
  period   <- generator$period
  cos1 <- cos(2 * pi * day_grid / period); sen1 <- sin(2 * pi * day_grid / period)
  cos2 <- cos(4 * pi * day_grid / period); sen2 <- sin(4 * pi * day_grid / period)

  sim_final <- matrix(0, nrow = n_days, ncol = length(generator$vars))
  colnames(sim_final) <- generator$vars

  for (v in generator$vars) {
    mod       <- generator$seasonal_models[[v]]
    res_unpad <- sim_pad[, v] * mod$sd_res + mod$mu_res
    saz <- mod$coefs["(Intercept)"] +
      mod$coefs["cos1"] * cos1 + mod$coefs["sen1"] * sen1 +
      mod$coefs["cos2"] * cos2 + mod$coefs["sen2"] * sen2
    sim_final[, v] <- saz + res_unpad
  }

  df_sim <- as.data.frame(sim_final)

  lapply(seq_len(n_mc), function(i) {
    idx <- ((i - 1L) * n_obs + 1L):(i * n_obs)
    df_sim[idx, , drop = FALSE]
  })
}
