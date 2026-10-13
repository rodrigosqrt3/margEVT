# =============================================================================
# generator.R
# Stationary VAR(p) generator of annual covariate trajectories, used as the
# parametric trajectory law (approach B) in marginalize().
#
# Workflow:
#   gen <- fit_var_generator(fit, data)
#   mc  <- simulate_covariates(gen, n_mc = 3000L, seed = 2024L)
#   marginalize(fit, data, approaches = c("A", "B", "C"), mc_sample = mc)
#
# Base variables are deseasonalized by OLS on Fourier harmonics and
# standardized; a stable VAR(p) without intercept is fitted to the
# standardized anomalies by least squares, with p chosen by the Schwarz
# criterion on a common estimation sample. Optionally, a user `derive`
# function recomputes model covariates (lags, rolling means, differences,
# standardizations) from simulated base variables, so that deterministic
# relations among covariates are preserved in every trajectory.
# =============================================================================

# Lagged regressor matrix [V_{t-1}, ..., V_{t-p}] for rows (start:n).
.var_regressors <- function(Y, p, start) {
  n <- nrow(Y)
  do.call(cbind, lapply(seq_len(p), function(i)
    Y[(start - i):(n - i), , drop = FALSE]))
}

# Least-squares VAR(p) without intercept on rows (p+1):n. Returns the K x Kp
# coefficient matrix and the residual matrix, or NULL if rank deficient.
.fit_var_ols <- function(Y, p) {
  n <- nrow(Y)
  if (n - p <= ncol(Y) * p) return(NULL)
  X <- .var_regressors(Y, p, p + 1L)
  R <- Y[(p + 1L):n, , drop = FALSE]
  qx <- qr(X)
  if (qx$rank < ncol(X)) return(NULL)
  B <- qr.coef(qx, R)
  list(A = t(B), residuals = R - X %*% B)
}

# Schwarz (SC) lag selection on a common sample, as in vars::VARselect():
# SC(p) = log det(Sigma_p) + log(T) / T * p * K^2, with T = n - lag_max.
.select_var_order <- function(resid_sel, lag_max) {
  Y <- as.matrix(resid_sel)
  n <- nrow(Y)
  K <- ncol(Y)
  T_eff <- n - lag_max
  R <- Y[(lag_max + 1L):n, , drop = FALSE]
  sc <- vapply(seq_len(lag_max), function(p) {
    X <- .var_regressors(Y, p, lag_max + 1L)
    qx <- qr(X)
    if (qx$rank < ncol(X)) return(NA_real_)
    E <- R - X %*% qr.coef(qx, R)
    d <- determinant(crossprod(E) / T_eff, logarithm = TRUE)
    as.numeric(d$modulus) + log(T_eff) / T_eff * p * K^2
  }, numeric(1L))
  if (all(is.na(sc))) return(NA_real_)
  which.min(sc)
}

.companion_matrix <- function(A, K, p) {
  if (p == 1L) return(A)
  lower <- cbind(diag(K * (p - 1L)), matrix(0, K * (p - 1L), K))
  rbind(A, lower)
}

# Moduli of the companion-matrix eigenvalues (inverse characteristic roots).
.var_root_moduli <- function(A, K, p) {
  Mod(eigen(.companion_matrix(A, K, p), only.values = TRUE)$values)
}

.stationary_state_covariance <- function(A, Sigma, K, p) {
  Fmat <- .companion_matrix(A, K, p)
  d <- K * p
  Q <- matrix(0, d, d)
  Q[seq_len(K), seq_len(K)] <- Sigma
  lhs <- diag(d * d) - kronecker(Fmat, Fmat)
  vec_p <- tryCatch(solve(lhs, as.vector(Q)), error = function(e) NULL)
  if (is.null(vec_p)) return(NULL)
  P <- matrix(vec_p, d, d)
  P <- (P + t(P)) / 2
  if (any(!is.finite(P))) NULL else P
}

# n independent draws from N(0, Sigma) for a PSD Sigma, as columns.
.rmvnorm_psd_cols <- function(Sigma, n) {
  eig <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow(Sigma))
  root %*% matrix(stats::rnorm(nrow(Sigma) * n), nrow(Sigma), n)
}

.regularize_covariance <- function(Sigma, relative_floor = 1e-10) {
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  scale_ref <- max(1, max(abs(eig$values)))
  values <- pmax(eig$values, relative_floor * scale_ref)
  out <- eig$vectors %*% (values * t(eig$vectors))
  (out + t(out)) / 2
}

#' Fit a stationary VAR generator of covariate trajectories
#'
#' Deseasonalizes each base variable by OLS on \code{n_harmonics} Fourier
#' harmonic pairs, standardizes the residuals, and fits a stable
#' vector autoregression VAR(p) without intercept by least squares, with the
#' order selected by the Schwarz criterion on a common estimation sample
#' (the convention of \code{vars::VARselect}). A univariate base variable
#' gives a stable AR(p) model.
#'
#' By default the base variables are the active covariates of \code{fit}.
#' When model covariates are deterministic transformations of a smaller set
#' of base variables (lags, rolling means, differences, interactions,
#' standardizations), supply the base variables in \code{vars} and a
#' \code{derive} function: \code{\link{simulate_covariates}} then simulates
#' the base variables and recomputes the model covariates on each
#' trajectory, preserving their deterministic relations.
#'
#' Seasonal harmonics present in \code{data} (\code{cos<k>}, \code{sen<k>} or
#' \code{sin<k>}) are reused; otherwise they are computed from the row index,
#' which assumes contiguous rows starting at the first position of a year.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param data Data frame containing the base variables (and, optionally,
#'   seasonal harmonic columns).
#' @param vars Character vector of base variables. Default \code{NULL}: the
#'   active covariates of \code{fit} present in \code{data}.
#' @param period Numeric seasonal period in observations. If \code{NULL},
#'   inherits \code{fit$obs_per_year}; legacy objects fall back to \code{365.25}.
#' @param lag_max Integer. Maximum VAR order considered. Default \code{NULL}:
#'   \code{min(5, floor(n / (3 K)))}.
#' @param n_harmonics Positive integer number of harmonic pairs used for
#'   deseasonalization. Default \code{2}.
#' @param derive Optional function taking a data frame of simulated base
#'   variables (plus the harmonic columns) for one trajectory and returning a
#'   data frame with the same number of rows containing the model covariates.
#' @param derive_lead Non-negative integer number of leading observations
#'   simulated before each year so that \code{derive} has the history it
#'   needs (e.g. the window of a lag or rolling mean). Default \code{0}.
#'
#' @return An object of class \code{nhpp_var_generator}, including the
#'   coefficient matrix, innovation covariance, root moduli and spectral
#'   radius, or \code{NULL} (with a warning) if there are no variables to
#'   model. An unstable selected VAR causes an error.
#'
#' @examples
#' set.seed(1)
#' n <- 3650
#' t <- rep(1:365, 10)
#' x <- as.numeric(stats::arima.sim(list(ar = 0.7), n)) + 2 * cos(2 * pi * t / 365.25)
#' df <- cbind(data.frame(x = x), seasonal_harmonics(t))
#' df$y <- 10 + df$x + rexp(n, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 20, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' gen <- fit_var_generator(fit, df)
#' gen$p_opt
#' gen$spectral_radius
#'
#' # Model covariates derived from a base variable: the lagged value of x
#' df$x_lag1 <- c(NA, head(df$x, -1))
#' gen_d <- fit_var_generator(
#'   fit, df[-1, ], vars = "x", derive_lead = 1,
#'   derive = function(b) data.frame(x = b$x, x_lag1 = c(NA, head(b$x, -1))))
#' head(simulate_covariates(gen_d, n_mc = 2, seed = 1)[[1]])
#'
#' @export
fit_var_generator <- function(fit, data, vars = NULL, period = NULL,
                              lag_max = NULL, n_harmonics = 2L,
                              derive = NULL, derive_lead = 0L) {

  if (!inherits(fit, "nhpp_fit"))
    stop("fit_var_generator: `fit` must be an nhpp_fit object.")
  if (!is.data.frame(data))
    stop("fit_var_generator: `data` must be a data frame.")

  period <- if (is.null(period)) {
    if (is.null(fit$obs_per_year)) 365.25 else fit$obs_per_year
  } else period
  if (!is.numeric(period) || length(period) != 1L ||
      !is.finite(period) || period <= 0)
    stop("fit_var_generator: `period` must be a single positive numeric value.")
  if (!is.null(lag_max) &&
      (!is.numeric(lag_max) || length(lag_max) != 1L ||
       !is.finite(lag_max) || lag_max < 1))
    stop("fit_var_generator: `lag_max` must be a positive integer.")
  if (!is.numeric(n_harmonics) || length(n_harmonics) != 1L ||
      !is.finite(n_harmonics) || n_harmonics < 1)
    stop("fit_var_generator: `n_harmonics` must be a positive integer.")
  n_harmonics <- as.integer(n_harmonics)
  if (!is.null(derive) && !is.function(derive))
    stop("fit_var_generator: `derive` must be a function or NULL.")
  if (!is.numeric(derive_lead) || length(derive_lead) != 1L ||
      !is.finite(derive_lead) || derive_lead < 0)
    stop("fit_var_generator: `derive_lead` must be a non-negative integer.")

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

  # Seasonal harmonic columns: reuse those in data, otherwise build them.
  sine <- if (all(paste0("sin", seq_len(n_harmonics)) %in% names(data)))
    "sin" else "sen"
  harm_names <- names(seasonal_harmonics(1, period, n_harmonics, sine))
  if (!all(harm_names %in% names(data))) {
    harm <- seasonal_harmonics(seq_len(nrow(data)), period, n_harmonics, sine)
    data[harm_names] <- harm
  }

  df_base <- data[, unique(c(vars, harm_names)), drop = FALSE]
  df_base <- df_base[stats::complete.cases(df_base), , drop = FALSE]
  if (nrow(df_base) < 10L)
    stop("fit_var_generator: too few complete rows after removing NAs.")

  seasonal_models <- list()
  stochastic <- logical(length(vars))
  resid_mat <- matrix(NA_real_, nrow = nrow(df_base), ncol = length(vars),
                      dimnames = list(NULL, vars))
  Xh <- cbind("(Intercept)" = 1, as.matrix(df_base[, harm_names, drop = FALSE]))
  for (i in seq_along(vars)) {
    v      <- vars[i]
    f      <- stats::lm.fit(Xh, df_base[[v]])
    res    <- as.numeric(f$residuals)
    mu_res <- mean(res)
    sd_res <- stats::sd(res)
    # A variable fully explained by its seasonal cycle is deterministic: it
    # is reconstructed from the seasonal fit and excluded from the VAR.
    stochastic[i] <- is.finite(sd_res) && sd_res >= 1e-8
    if (!stochastic[i]) sd_res <- 1
    coefs <- f$coefficients
    coefs[is.na(coefs)] <- 0
    seasonal_models[[v]] <- list(coefs = stats::setNames(coefs, colnames(Xh)),
                                 mu_res = mu_res, sd_res = sd_res)
    resid_mat[, v] <- (res - mu_res) / sd_res
  }

  stoch_vars <- vars[stochastic]
  K <- length(stoch_vars)
  if (K == 0L) {
    p_opt <- 0L
    A <- matrix(0, 0L, 0L)
    residual_matrix <- matrix(0, 0L, 0L)
    Sigma <- matrix(0, 0L, 0L)
    root_moduli <- numeric(0L)
    spectral_radius <- 0
  } else {
    Y <- resid_mat[, stoch_vars, drop = FALSE]
    if (is.null(lag_max))
      lag_max <- max(1L, min(5L, floor(nrow(Y) / (3L * K))))
    lag_max <- min(as.integer(lag_max), max(1L, nrow(Y) - 2L))

    p_opt <- as.integer(.select_var_order(Y, lag_max))
    if (!is.finite(p_opt))
      stop("fit_var_generator: VAR lag selection did not return a finite order.")
    ols <- .fit_var_ols(Y, p_opt)
    if (is.null(ols))
      stop("fit_var_generator: the selected VAR is not identifiable (rank deficient).")
    A <- ols$A
    residual_matrix <- ols$residuals
    Sigma <- crossprod(residual_matrix) / nrow(residual_matrix)
    root_moduli <- .var_root_moduli(A, K, p_opt)
    spectral_radius <- max(root_moduli)
    if (!is.finite(spectral_radius) || spectral_radius >= 1)
      stop(sprintf(
        paste0(
          "fit_var_generator: selected VAR(%d) is unstable ",
          "(spectral radius = %.6f). Reduce `lag_max`, revise the covariate set, ",
          "or model the trajectory law separately."
        ),
        p_opt, spectral_radius
      ))
  }

  structure(
    list(
      vars            = vars,
      stochastic_vars = stoch_vars,
      seasonal_models = seasonal_models,
      harmonics       = harm_names,
      n_harmonics     = n_harmonics,
      A               = A,
      Sigma           = Sigma,
      residuals       = residual_matrix,
      K               = K,
      p_opt           = p_opt,
      root_moduli     = root_moduli,
      spectral_radius = spectral_radius,
      period          = period,
      var_colnames    = vars,
      derive          = derive,
      derive_lead     = as.integer(derive_lead)
    ),
    class = "nhpp_var_generator"
  )
}


# Simulates n_mc independent VAR(p) anomaly paths of length n_total at once.
# Returns an array [n_total, K, n_mc].
.simulate_var_paths <- function(A, chol_S, stationary_covariance, K, p, n_mc,
                                n_total, initialization, burn_in, keep_lag,
                                seed) {
  .maybe_seeded(seed, {
    state <- if (initialization == "stationary")
      .rmvnorm_psd_cols(stationary_covariance, n_mc) else
        matrix(0, K * p, n_mc)
    step <- function(state) {
      new_y <- A %*% state + chol_S %*% matrix(stats::rnorm(K * n_mc), K, n_mc)
      list(y = new_y, state = rbind(new_y, state[keep_lag, , drop = FALSE]))
    }
    if (initialization == "burnin")
      for (b in seq_len(burn_in)) state <- step(state)$state
    out <- array(NA_real_, dim = c(n_total, K, n_mc))
    for (tt in seq_len(n_total)) {
      s <- step(state)
      out[tt, , ] <- s$y
      state <- s$state
    }
    out
  })
}

#' Simulate stationary covariate trajectories from a VAR generator
#'
#' Simulates \code{n_mc} independent annual trajectories. All trajectories
#' are propagated simultaneously, so the cost grows with the trajectory
#' length rather than with \code{n_mc} times the length. Each trajectory is
#' initialized from the stationary Gaussian law of the companion state
#' (default) or from zero followed by \code{burn_in} transitions. Seasonality
#' is re-added at positions \code{1, ..., n_obs}; if the generator has a
#' \code{derive} function, model covariates are recomputed from the
#' simulated base variables, using \code{derive_lead} additional leading
#' observations that are discarded afterwards.
#'
#' @param generator Object from \code{\link{fit_var_generator}}. If
#'   \code{NULL}, returns \code{NULL} (stationary/no-covariate model).
#' @param n_mc Integer. Number of independently simulated annual trajectories.
#' @param n_obs Integer observations per simulated year. If \code{NULL},
#'   uses the rounded seasonal period stored in \code{generator}.
#' @param burn_in Integer. Burn-in applied when
#'   \code{initialization = "burnin"}. Default \code{300L}.
#' @param initialization Character. \code{"stationary"} (default) draws the
#'   initial companion state from the fitted stationary Gaussian law;
#'   \code{"burnin"} starts at zero and applies \code{burn_in} transitions.
#' @param seed Integer seed or \code{NULL}. When supplied, the caller's
#'   random-number stream is preserved. Default \code{NULL} (use the current
#'   stream).
#'
#' @return A list of length \code{n_mc}, each element an \code{n_obs}-row
#'   data frame of covariates, directly usable as \code{mc_sample} in
#'   \code{\link{marginalize}} (approach B).
#'
#' @examples
#' set.seed(1)
#' n <- 3650
#' t <- rep(1:365, 10)
#' df <- cbind(data.frame(x = as.numeric(stats::arima.sim(list(ar = 0.6), n))),
#'             seasonal_harmonics(t))
#' df$y <- 10 + df$x + rexp(n, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 20, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' mc <- simulate_covariates(fit_var_generator(fit, df), n_mc = 3, seed = 1)
#' length(mc); dim(mc[[1]])
#'
#' @export
simulate_covariates <- function(generator, n_mc, n_obs = NULL,
                                burn_in = 300L,
                                initialization = c("stationary", "burnin"),
                                seed = NULL) {

  if (is.null(generator)) return(NULL)
  if (!inherits(generator, "nhpp_var_generator"))
    stop("simulate_covariates: `generator` must come from fit_var_generator().")
  if (!is.numeric(n_mc) || length(n_mc) != 1L ||
      !is.finite(n_mc) || n_mc < 1)
    stop("simulate_covariates: `n_mc` must be a positive integer.")
  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
      !is.finite(burn_in) || burn_in < 0)
    stop("simulate_covariates: `burn_in` must be a non-negative integer.")
  .validate_seed(seed, "simulate_covariates")
  n_mc    <- as.integer(n_mc)
  burn_in <- as.integer(burn_in)
  initialization <- match.arg(initialization)

  inherited_period <- if (is.null(generator$period)) 365.25 else generator$period
  n_obs  <- if (is.null(n_obs)) as.integer(round(inherited_period)) else
    as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L)
    stop("simulate_covariates: `n_obs` must be a positive integer.")

  vars  <- generator$vars
  stoch <- if (is.null(generator$stochastic_vars)) vars else
    generator$stochastic_vars
  K <- generator$K
  p <- generator$p_opt
  A <- generator$A
  lead    <- if (is.null(generator$derive_lead)) 0L else generator$derive_lead
  n_total <- n_obs + lead

  if (K == 0L) {
    anom <- array(0, dim = c(n_total, 0L, n_mc))
  } else {
    Sigma_raw <- if (is.null(generator$Sigma))
      crossprod(generator$residuals) / nrow(generator$residuals) else
        generator$Sigma
    Sigma  <- .regularize_covariance(Sigma_raw)
    chol_S <- t(chol(Sigma))
    stationary_covariance <- if (initialization == "stationary")
      .stationary_state_covariance(A, Sigma, K, p) else NULL
    if (initialization == "stationary" && is.null(stationary_covariance))
      stop("simulate_covariates: stationary initialization failed; use initialization='burnin' for a sensitivity analysis.")
    keep_lag <- seq_len(K * (p - 1L))
    anom <- .simulate_var_paths(A, chol_S, stationary_covariance, K, p,
                                n_mc, n_total, initialization, burn_in,
                                keep_lag, seed)
  }

  # Seasonal re-injection at positions (1 - lead), ..., n_obs.
  t_pos <- seq_len(n_total) - lead
  harm_names <- if (is.null(generator$harmonics))
    c("cos1", "sen1", "cos2", "sen2") else generator$harmonics
  Xh <- cbind("(Intercept)" = 1,
              as.matrix(as.data.frame(stats::setNames(
                lapply(harm_names, .harmonic_column, t = t_pos,
                       period = inherited_period), harm_names))))

  seasonal <- vapply(vars, function(v) {
    mod <- generator$seasonal_models[[v]]
    as.numeric(Xh %*% mod$coefs[colnames(Xh)]) + mod$mu_res
  }, numeric(n_total))
  seasonal <- matrix(seasonal, nrow = n_total)
  sd_res <- vapply(vars, function(v) generator$seasonal_models[[v]]$sd_res,
                   numeric(1L))
  stoch_idx <- match(stoch, vars)
  base_paths <- lapply(seq_len(n_mc), function(r) {
    vals <- seasonal
    if (K > 0L)
      vals[, stoch_idx] <- vals[, stoch_idx] +
        sweep(matrix(anom[, , r], nrow = n_total), 2L, sd_res[stoch_idx], `*`)
    colnames(vals) <- vars
    data.frame(vals, check.names = FALSE)
  })

  if (is.null(generator$derive)) return(base_paths)

  keep <- (lead + 1L):n_total
  lapply(base_paths, function(b) {
    b_full <- cbind(b, as.data.frame(Xh[, -1L, drop = FALSE]))
    d <- generator$derive(b_full)
    if (!is.data.frame(d) || nrow(d) != n_total)
      stop("simulate_covariates: `derive` must return a data frame with ",
           n_total, " rows.")
    d <- d[keep, , drop = FALSE]
    if (any(!vapply(d, function(col) is.numeric(col) && all(is.finite(col)),
                    logical(1L))))
      stop("simulate_covariates: `derive` returned non-numeric or non-finite ",
           "values after the leading observations; increase `derive_lead`.")
    rownames(d) <- NULL
    d
  })
}
