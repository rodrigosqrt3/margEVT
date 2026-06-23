# =============================================================================
# likelihood.R
# Penalized negative log-likelihood and exact analytical gradient for the
# non-homogeneous point process extreme value model.
# =============================================================================

#' Penalized negative log-likelihood for the NHPP model
#'
#' @param par Numeric vector of parameters: c(beta_mu, beta_sigma, beta_xi).
#' @param dm Design matrix list from \code{\link{build_design_matrices}}.
#' @param y Numeric vector of observations.
#' @param threshold Numeric scalar. The extreme value threshold u.
#' @param lambda Penalty strength. Either a scalar or named vector
#'   \code{c(mu=, sigma=, xi=)}.
#' @param alpha Elastic-net mixing parameter in [0,1]. 1 = LASSO, 0 = ridge.
#'   Either a scalar or named vector \code{c(mu=, sigma=, xi=)}.
#' @param pen_xi Logical. Should the shape parameter covariates be penalized?
#' @param eps_l1 Small positive constant for smooth L1 approximation.
#' @param obs_per_year Number of observations per year (e.g. 365.25 for daily).
#'
#' @return Scalar. The penalized negative log-likelihood value.
#'
#' @export
pp_nllh <- function(par, dm, y, threshold,
                    lambda = 0.1, alpha = 0.5,
                    pen_xi = TRUE, eps_l1 = 1e-10,
                    obs_per_year = 365.25) {

  X_mu  <- dm$X_mu
  X_sig <- dm$X_sigma
  X_xi  <- dm$X_xi
  p_mu  <- ncol(X_mu)
  p_sig <- ncol(X_sig)
  p_xi  <- ncol(X_xi)

  beta_mu  <- par[seq_len(p_mu)]
  beta_sig <- par[p_mu + seq_len(p_sig)]
  beta_xi  <- par[p_mu + p_sig + seq_len(p_xi)]

  mu_t  <- as.numeric(X_mu  %*% beta_mu)
  sig_t <- exp(as.numeric(X_sig %*% beta_sig))
  xi_t  <- as.numeric(X_xi  %*% beta_xi)

  z_u <- 1 + xi_t * (threshold - mu_t) / sig_t
  if (any(z_u <= 0, na.rm = TRUE)) return(1e9)

  gumbel <- abs(xi_t) < 1e-6
  tau    <- numeric(length(mu_t))
  tau[ gumbel] <- exp(pmax(-500, -(threshold - mu_t[gumbel]) / sig_t[gumbel]))
  tau[!gumbel] <- z_u[!gumbel]^(-1 / xi_t[!gumbel])

  Lambda <- sum(tau) / obs_per_year

  exc_idx <- which(y > threshold)
  log_f   <- 0
  if (length(exc_idx) > 0L) {
    y_e   <- y[exc_idx]
    mu_e  <- mu_t[exc_idx]
    sig_e <- sig_t[exc_idx]
    xi_e  <- xi_t[exc_idx]
    z_e   <- 1 + xi_e * (y_e - mu_e) / sig_e
    if (any(z_e <= 0, na.rm = TRUE)) return(1e9)
    gum_e <- abs(xi_e) < 1e-6
    lf    <- numeric(length(y_e))
    lf[ gum_e] <- -log(sig_e[gum_e])  - (y_e[gum_e]  - mu_e[gum_e])  / sig_e[gum_e]
    lf[!gum_e] <- -log(sig_e[!gum_e]) - (1 + 1 / xi_e[!gum_e]) * log(z_e[!gum_e])
    log_f <- sum(lf)
  }

  nllh_raw <- Lambda - log_f

  lam <- function(k) if (length(lambda) > 1L) lambda[k] else lambda
  alp <- function(k) if (length(alpha)  > 1L) alpha[k]  else alpha

  pen <- 0
  if (length(dm$idx_pen_mu) > 0L) {
    g <- beta_mu[dm$idx_pen_mu]
    pen <- pen + lam("mu") * sum(alp("mu") * sqrt(g^2 + eps_l1) +
                                   0.5 * (1 - alp("mu")) * g^2)
  }
  if (length(dm$idx_pen_sigma) > 0L) {
    g <- beta_sig[dm$idx_pen_sigma]
    pen <- pen + lam("sigma") * sum(alp("sigma") * sqrt(g^2 + eps_l1) +
                                      0.5 * (1 - alp("sigma")) * g^2)
  }
  if (pen_xi && length(dm$idx_pen_xi) > 0L) {
    g <- beta_xi[dm$idx_pen_xi]
    pen <- pen + lam("xi") * sum(alp("xi") * sqrt(g^2 + eps_l1) +
                                   0.5 * (1 - alp("xi")) * g^2)
  }

  nllh_raw + pen
}


#' Analytical gradient of the penalized negative log-likelihood
#'
#' @inheritParams pp_nllh
#'
#' @return Numeric vector of length \code{length(par)}.
#'
#' @export
pp_grad <- function(par, dm, y, threshold,
                    lambda = 0.1, alpha = 0.5,
                    pen_xi = TRUE, eps_l1 = 1e-10,
                    obs_per_year = 365.25) {

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
  dL_dmu[valid &  gumbel] <- tau[valid &  gumbel] / sig_t[valid &  gumbel]
  dL_dmu[valid & !gumbel] <- z_u[valid & !gumbel]^(
    -1 / xi_t[valid & !gumbel] - 1) / sig_t[valid & !gumbel]
  g_mu <- g_mu + dL_dmu / obs_per_year

  wu <- (threshold - mu_t) / sig_t
  dL_dlsig <- numeric(n)
  dL_dlsig[valid &  gumbel] <- tau[valid &  gumbel] * wu[valid &  gumbel]
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

  grad_mu    <- as.numeric(t(X_mu)    %*% g_mu)
  grad_sigma <- as.numeric(t(X_sigma) %*% g_lsig)
  grad_xi    <- as.numeric(t(X_xi)    %*% g_xi)

  lam <- function(k) if (length(lambda) > 1L) lambda[k] else lambda
  alp <- function(k) if (length(alpha)  > 1L) alpha[k]  else alpha

  pen_grad <- function(gamma, lam_k, alp_k) {
    if (lam_k == 0 || length(gamma) == 0L) return(rep(0, length(gamma)))
    lam_k * (alp_k * gamma / sqrt(gamma^2 + eps_l1) + (1 - alp_k) * gamma)
  }

  if (length(dm$idx_pen_mu) > 0L)
    grad_mu[dm$idx_pen_mu] <- grad_mu[dm$idx_pen_mu] +
    pen_grad(beta_mu[dm$idx_pen_mu], lam("mu"), alp("mu"))

  if (length(dm$idx_pen_sigma) > 0L)
    grad_sigma[dm$idx_pen_sigma] <- grad_sigma[dm$idx_pen_sigma] +
    pen_grad(beta_sigma[dm$idx_pen_sigma], lam("sigma"), alp("sigma"))

  if (pen_xi && length(dm$idx_pen_xi) > 0L)
    grad_xi[dm$idx_pen_xi] <- grad_xi[dm$idx_pen_xi] +
    pen_grad(beta_xi[dm$idx_pen_xi], lam("xi"), alp("xi"))

  c(grad_mu, grad_sigma, grad_xi)
}
