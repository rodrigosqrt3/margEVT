# =============================================================================
# likelihood.R
# Penalized negative log-likelihood and exact analytical gradient for the
# non-homogeneous point process extreme value model.
#
# Contents:
#   - .resolve_block()        : block-level control lookup (internal)
#   - .penalty_vectors()      : coordinate-wise lambda and alpha (internal)
#   - .nhpp_eval()            : unpenalized value and gradient (internal)
#   - .penalty_value_grad()   : smooth elastic-net value and gradient (internal)
#   - .objective_factory()    : cached fn/gr pair for optim() (internal)
#   - .kkt_violation()        : KKT-based stationarity measure (internal)
#   - pp_nllh(), pp_grad()    : exported value and gradient
# =============================================================================

# Block-level control lookup. Accepts a scalar, a named vector c(mu=, sigma=,
# xi=), or an unnamed length-3 vector in the order mu, sigma, xi.
.resolve_block <- function(x, block) {
  if (length(x) == 1L) return(unname(x))
  if (!is.null(names(x)) && block %in% names(x)) return(unname(x[[block]]))
  unname(x[[match(block, c("mu", "sigma", "xi"))]])
}

# Coordinate-wise penalty strength and mixing parameter. Unpenalized
# coordinates receive lambda = 0. Optional penalty weights multiply lambda
# coordinate by coordinate (adaptive LASSO and user-supplied weights).
.penalty_vectors <- function(dm, lambda, alpha, pen_xi,
                             penalty_weights = NULL) {
  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)
  p     <- p_mu + p_sig + p_xi
  lam <- numeric(p)
  alp <- numeric(p)

  idx_mu  <- dm$idx_pen_mu
  idx_sig <- p_mu + dm$idx_pen_sigma
  idx_xi  <- if (isTRUE(pen_xi)) p_mu + p_sig + dm$idx_pen_xi else integer(0L)

  lam[idx_mu]  <- .resolve_block(lambda, "mu")
  lam[idx_sig] <- .resolve_block(lambda, "sigma")
  lam[idx_xi]  <- .resolve_block(lambda, "xi")
  alp[idx_mu]  <- .resolve_block(alpha, "mu")
  alp[idx_sig] <- .resolve_block(alpha, "sigma")
  alp[idx_xi]  <- .resolve_block(alpha, "xi")

  if (!is.null(penalty_weights)) lam <- lam * as.numeric(penalty_weights)
  list(lam = lam, alp = alp)
}

# Unpenalized negative log-likelihood and its exact gradient.
# Returns value = 1e9 outside the point-process support; the gradient is
# still evaluated with the support masks, as optimizers only use gradients
# at accepted (finite-objective) iterates.
.nhpp_eval <- function(par, dm, y, threshold, obs_per_year = 365.25,
                       grad = TRUE) {
  X_mu    <- dm$X_mu
  X_sigma <- dm$X_sigma
  X_xi    <- dm$X_xi
  p_mu    <- ncol(X_mu)
  p_sig   <- ncol(X_sigma)
  p_xi    <- ncol(X_xi)
  n       <- nrow(X_mu)

  beta_mu    <- par[seq_len(p_mu)]
  beta_sigma <- par[p_mu + seq_len(p_sig)]
  beta_xi    <- par[p_mu + p_sig + seq_len(p_xi)]

  mu_t  <- as.numeric(X_mu    %*% beta_mu)
  sig_t <- exp(as.numeric(X_sigma %*% beta_sigma))
  xi_t  <- as.numeric(X_xi    %*% beta_xi)

  xi_tol  <- 1e-6
  exc_idx <- which(y > threshold)

  # ---- value --------------------------------------------------------------
  value <- NA_real_
  if (any(!is.finite(mu_t)) || any(!is.finite(sig_t)) ||
      any(sig_t <= 0) || any(!is.finite(xi_t))) {
    value <- 1e9
  } else {
    z_u <- 1 + xi_t * (threshold - mu_t) / sig_t
    if (any(!is.finite(z_u)) || any(z_u <= 0)) {
      value <- 1e9
    } else {
      gumbel <- abs(xi_t) < xi_tol
      tau    <- numeric(n)
      tau[ gumbel] <- exp(pmax(-500, -(threshold - mu_t[gumbel]) / sig_t[gumbel]))
      tau[!gumbel] <- z_u[!gumbel]^(-1 / xi_t[!gumbel])
      Lambda <- sum(tau) / obs_per_year
      log_f  <- 0
      if (length(exc_idx) > 0L) {
        y_e   <- y[exc_idx]
        mu_e  <- mu_t[exc_idx]
        sig_e <- sig_t[exc_idx]
        xi_e  <- xi_t[exc_idx]
        z_e   <- 1 + xi_e * (y_e - mu_e) / sig_e
        if (any(!is.finite(z_e)) || any(z_e <= 0)) {
          value <- 1e9
        } else {
          gum_e <- abs(xi_e) < xi_tol
          lf    <- numeric(length(y_e))
          lf[ gum_e] <- -log(sig_e[gum_e]) - (y_e[gum_e] - mu_e[gum_e]) / sig_e[gum_e]
          lf[!gum_e] <- -log(sig_e[!gum_e]) - (1 + 1 / xi_e[!gum_e]) * log(z_e[!gum_e])
          log_f <- sum(lf)
        }
      }
      if (is.na(value)) value <- Lambda - log_f
    }
  }
  if (!isTRUE(grad)) return(list(value = value, grad = NULL))

  # ---- gradient -----------------------------------------------------------
  g_mu   <- numeric(n)
  g_lsig <- numeric(n)
  g_xi   <- numeric(n)

  z_u    <- 1 + xi_t * (threshold - mu_t) / sig_t
  valid  <- is.finite(z_u) & (z_u > 0)
  gumbel <- abs(xi_t) < xi_tol

  tau <- numeric(n)
  tau[valid &  gumbel] <- exp(pmax(-500,
    -(threshold - mu_t[valid & gumbel]) / sig_t[valid & gumbel]))
  tau[valid & !gumbel] <- z_u[valid & !gumbel]^(-1 / xi_t[valid & !gumbel])

  dL_dmu <- numeric(n)
  dL_dmu[valid &  gumbel] <- tau[valid & gumbel] / sig_t[valid & gumbel]
  dL_dmu[valid & !gumbel] <- z_u[valid & !gumbel]^(
    -1 / xi_t[valid & !gumbel] - 1) / sig_t[valid & !gumbel]
  g_mu <- g_mu + dL_dmu / obs_per_year

  wu <- (threshold - mu_t) / sig_t
  dL_dlsig <- numeric(n)
  dL_dlsig[valid &  gumbel] <- tau[valid & gumbel] * wu[valid & gumbel]
  dL_dlsig[valid & !gumbel] <- z_u[valid & !gumbel]^(
    -1 / xi_t[valid & !gumbel] - 1) * wu[valid & !gumbel]
  g_lsig <- g_lsig + dL_dlsig / obs_per_year

  dL_dxi <- numeric(n)
  dL_dxi[valid & gumbel] <- 0.5 * tau[valid & gumbel] * wu[valid & gumbel]^2
  ok_xi <- valid & !gumbel
  if (any(ok_xi)) {
    dL_dxi[ok_xi] <- tau[ok_xi] * (
      log(z_u[ok_xi]) / xi_t[ok_xi]^2 -
        wu[ok_xi] / (z_u[ok_xi] * xi_t[ok_xi])
    )
  }
  g_xi <- g_xi + dL_dxi / obs_per_year

  if (length(exc_idx) > 0L) {
    y_e   <- y[exc_idx]
    mu_e  <- mu_t[exc_idx]
    sig_e <- sig_t[exc_idx]
    xi_e  <- xi_t[exc_idx]
    z_e   <- pmax(1 + xi_e * (y_e - mu_e) / sig_e, 1e-6)
    gum_e <- abs(xi_e) < xi_tol
    we    <- (y_e - mu_e) / sig_e

    dn_dmu <- numeric(length(exc_idx))
    dn_dmu[ gum_e] <- -1 / sig_e[gum_e]
    dn_dmu[!gum_e] <- -(1 + xi_e[!gum_e]) / (sig_e[!gum_e] * z_e[!gum_e])
    g_mu[exc_idx] <- g_mu[exc_idx] + dn_dmu

    dn_dlsig <- numeric(length(exc_idx))
    dn_dlsig[ gum_e] <- 1 - we[gum_e]
    dn_dlsig[!gum_e] <- 1 - (1 + xi_e[!gum_e]) * we[!gum_e] / z_e[!gum_e]
    g_lsig[exc_idx] <- g_lsig[exc_idx] + dn_dlsig

    dn_dxi <- numeric(length(exc_idx))
    dn_dxi[gum_e] <- we[gum_e] - 0.5 * we[gum_e]^2
    if (any(!gum_e)) {
      dn_dxi[!gum_e] <- -log(z_e[!gum_e]) / xi_e[!gum_e]^2 +
        (1 + 1 / xi_e[!gum_e]) * we[!gum_e] / z_e[!gum_e]
    }
    g_xi[exc_idx] <- g_xi[exc_idx] + dn_dxi
  }

  list(value = value,
       grad  = c(as.numeric(crossprod(X_mu,    g_mu)),
                 as.numeric(crossprod(X_sigma, g_lsig)),
                 as.numeric(crossprod(X_xi,    g_xi))))
}

# Smooth elastic-net penalty value and gradient for coordinate vectors.
.penalty_value_grad <- function(par, pv, eps_l1 = 1e-10) {
  k <- which(pv$lam > 0)
  grad <- numeric(length(par))
  if (length(k) == 0L) return(list(value = 0, grad = grad))
  g <- par[k]
  value <- sum(pv$lam[k] * (pv$alp[k] * sqrt(g^2 + eps_l1) +
                              0.5 * (1 - pv$alp[k]) * g^2))
  grad[k] <- pv$lam[k] * (pv$alp[k] * g / sqrt(g^2 + eps_l1) +
                            (1 - pv$alp[k]) * g)
  list(value = value, grad = grad)
}

# Cached objective for optim(): value and gradient are computed together and
# reused when optim() requests both at the same parameter vector.
.objective_factory <- function(dm, y, threshold, pv, eps_l1 = 1e-10,
                               obs_per_year = 365.25) {
  last_par <- NULL
  last     <- NULL
  evaluate <- function(par) {
    if (is.null(last_par) || !identical(par, last_par)) {
      core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = TRUE)
      pen  <- .penalty_value_grad(par, pv, eps_l1)
      value <- if (core$value >= 1e9) 1e9 else core$value + pen$value
      last <<- list(value = value, grad = core$grad + pen$grad)
      last_par <<- par
    }
    last
  }
  list(fn = function(par) evaluate(par)$value,
       gr = function(par) evaluate(par)$grad)
}

# KKT-based stationarity measure for the smoothed elastic-net objective.
# For penalized coordinates numerically at zero (|theta| < active_tol), the
# subgradient condition |d(-loglik + ridge)/d theta| <= lambda * alpha is
# checked; every other coordinate must have a near-zero gradient.
.kkt_violation <- function(par, grad_total, pv, active_tol, eps_l1 = 1e-10) {
  if (any(!is.finite(grad_total))) return(Inf)
  v <- abs(grad_total)
  at_zero <- pv$lam > 0 & pv$alp > 0 & abs(par) < active_tol
  if (any(at_zero)) {
    l1_part <- pv$lam * pv$alp * par / sqrt(par^2 + eps_l1)
    g_smooth <- grad_total - l1_part
    v[at_zero] <- pmax(0, abs(g_smooth[at_zero]) -
                         pv$lam[at_zero] * pv$alp[at_zero])
  }
  max(v)
}


#' Penalized negative log-likelihood for the NHPP model
#'
#' Evaluates the negative log-likelihood of the non-homogeneous Poisson
#' process (point-process) representation of threshold exceedances, plus a
#' smoothed elastic-net penalty on the penalized coefficients.
#'
#' @param par Numeric vector of parameters: c(beta_mu, beta_sigma, beta_xi).
#' @param dm Design matrix list from \code{\link{build_design_matrices}}.
#' @param y Numeric vector of observations.
#' @param threshold Numeric scalar. The extreme value threshold u.
#' @param lambda Penalty strength. Either a scalar or a vector
#'   \code{c(mu=, sigma=, xi=)} (named, or unnamed in that order).
#' @param alpha Elastic-net mixing parameter in [0,1]. 1 = LASSO, 0 = ridge.
#'   Either a scalar or a vector \code{c(mu=, sigma=, xi=)}.
#' @param pen_xi Logical. Should the shape parameter covariates be penalized?
#' @param eps_l1 Small positive constant for smooth L1 approximation.
#' @param obs_per_year Number of observations per year (e.g. 365.25 for daily).
#' @param penalty_weights Optional non-negative numeric vector of length
#'   \code{length(par)} multiplying \code{lambda} coordinate by coordinate
#'   (e.g. adaptive-LASSO weights). \code{NULL} (default) means unit weights.
#'
#' @return Scalar. The penalized negative log-likelihood value. Parameter
#'   vectors outside the point-process support return \code{1e9}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' dm <- build_design_matrices(df, loc_vars = "x")
#' par <- c(6, 0.1, log(1.5), 0.1)
#' pp_nllh(par, dm, df$y, threshold = 4, lambda = 0)
#' pp_nllh(par, dm, df$y, threshold = 4, lambda = 1)
#'
#' @export
pp_nllh <- function(par, dm, y, threshold,
                    lambda = 0.1, alpha = 0.5,
                    pen_xi = TRUE, eps_l1 = 1e-10,
                    obs_per_year = 365.25,
                    penalty_weights = NULL) {
  core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = FALSE)
  if (core$value >= 1e9) return(1e9)
  pv <- .penalty_vectors(dm, lambda, alpha, pen_xi, penalty_weights)
  core$value + .penalty_value_grad(par, pv, eps_l1)$value
}


#' Analytical gradient of the penalized negative log-likelihood
#'
#' @inheritParams pp_nllh
#'
#' @return Numeric vector of length \code{length(par)}.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' dm <- build_design_matrices(df, loc_vars = "x")
#' pp_grad(c(6, 0.1, log(1.5), 0.1), dm, df$y, threshold = 4, lambda = 0)
#'
#' @export
pp_grad <- function(par, dm, y, threshold,
                    lambda = 0.1, alpha = 0.5,
                    pen_xi = TRUE, eps_l1 = 1e-10,
                    obs_per_year = 365.25,
                    penalty_weights = NULL) {
  core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = TRUE)
  pv <- .penalty_vectors(dm, lambda, alpha, pen_xi, penalty_weights)
  core$grad + .penalty_value_grad(par, pv, eps_l1)$grad
}
