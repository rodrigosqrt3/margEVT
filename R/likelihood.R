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
#   - .shape_bound_map()      : box reparametrization of the shape bound (internal)
#   - .kkt_violation_bounded(): KKT measure under the shape bound (internal)
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
# The integrated intensity follows the point-process definition
# Lambda = sum_t [1 + xi_t (u - mu_t) / sigma_t]_+^(-1/xi_t) / n_y:
# when xi_t < 0 and the upper endpoint mu_t - sigma_t / xi_t lies below u,
# the intensity vanishes above u and period t contributes zero. The value is
# 1e9 only outside the model support: an exceedance beyond the upper endpoint,
# a threshold below the lower endpoint (xi_t > 0), or xi_t < shape_lower.
#
# This is the fast path, used for every evaluation. It computes z_u, log z_u
# and tau once and obtains z_u^(-1/xi - 1) as tau / z_u. Outside the support
# it returns 1e9 directly when only the value is requested (line-search
# trials); with grad = TRUE, and whenever some |xi_t| < 1e-6 (Gumbel limits),
# it defers to .nhpp_eval_ref(), so the masked gradient and the Gumbel
# expressions are unchanged. Both agree to rounding error (tests/test-v050.R).
.nhpp_eval <- function(par, dm, y, threshold, obs_per_year = 365.25,
                       grad = TRUE, shape_lower = -Inf) {
  X_mu <- dm$X_mu; X_sigma <- dm$X_sigma; X_xi <- dm$X_xi
  p_mu <- ncol(X_mu); p_sig <- ncol(X_sigma); p_xi <- ncol(X_xi)
  mu_t <- drop(X_mu %*% par[seq_len(p_mu)])
  lsig <- drop(X_sigma %*% par[p_mu + seq_len(p_sig)])
  sig_t <- exp(lsig)
  xi_t <- drop(X_xi %*% par[p_mu + p_sig + seq_len(p_xi)])

  reference <- function()
    .nhpp_eval_ref(par, dm, y, threshold, obs_per_year, grad, shape_lower)
  outside <- function()
    if (isTRUE(grad)) reference() else list(value = 1e9, grad = NULL)

  if (!all(is.finite(xi_t)) || any(abs(xi_t) < 1e-6)) return(reference())
  if (!all(is.finite(mu_t)) || !all(is.finite(sig_t)) || any(sig_t <= 0) ||
      any(xi_t < shape_lower - 1e-10))
    return(outside())

  wu  <- (threshold - mu_t) / sig_t
  z_u <- 1 + xi_t * wu
  if (!all(is.finite(z_u)) || any(xi_t > 0 & z_u <= 0)) return(outside())

  exc  <- which(y > threshold)
  mu_e <- mu_t[exc]; sig_e <- sig_t[exc]; xi_e <- xi_t[exc]
  we   <- (y[exc] - mu_e) / sig_e
  z_e  <- 1 + xi_e * we
  if (!all(is.finite(z_e)) || any(z_e <= 0)) return(outside())
  lze <- log(z_e)

  pos    <- z_u > 0
  allpos <- all(pos)
  if (allpos) {
    lz  <- log(z_u)
    tau <- exp(-lz / xi_t)
  } else {
    lz  <- numeric(length(z_u)); lz[pos]  <- log(z_u[pos])
    tau <- numeric(length(z_u)); tau[pos] <- exp(-lz[pos] / xi_t[pos])
  }
  value <- sum(tau) / obs_per_year - sum(-lsig[exc] - (1 + 1 / xi_e) * lze)
  if (!isTRUE(grad)) return(list(value = value, grad = NULL))

  if (allpos) {
    tz     <- tau / z_u
    g_mu   <- tz / sig_t
    g_lsig <- tz * wu
    g_xi   <- tau * (lz / xi_t^2 - wu / (z_u * xi_t))
  } else {
    g_mu <- g_lsig <- g_xi <- numeric(length(z_u))
    k  <- which(pos)
    tz <- tau[k] / z_u[k]
    g_mu[k]   <- tz / sig_t[k]
    g_lsig[k] <- tz * wu[k]
    g_xi[k]   <- tau[k] * (lz[k] / xi_t[k]^2 - wu[k] / (z_u[k] * xi_t[k]))
  }
  g_mu   <- g_mu / obs_per_year
  g_lsig <- g_lsig / obs_per_year
  g_xi   <- g_xi / obs_per_year
  a <- (1 + xi_e) / z_e
  g_mu[exc]   <- g_mu[exc]   - a / sig_e
  g_lsig[exc] <- g_lsig[exc] + 1 - a * we
  g_xi[exc]   <- g_xi[exc]   - lze / xi_e^2 + (1 + 1 / xi_e) * we / z_e

  list(value = value,
       grad  = c(as.numeric(crossprod(X_mu,    g_mu)),
                 as.numeric(crossprod(X_sigma, g_lsig)),
                 as.numeric(crossprod(X_xi,    g_xi))))
}

# Reference implementation (used by .nhpp_eval() for support violations with
# grad = TRUE and for Gumbel coordinates). The gradient uses the support
# masks; tau_t and its derivatives vanish continuously at z_u = 0 when
# -1 < xi_t < 0.
.nhpp_eval_ref <- function(par, dm, y, threshold, obs_per_year = 365.25,
                           grad = TRUE, shape_lower = -Inf) {
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
      any(sig_t <= 0) || any(!is.finite(xi_t)) ||
      any(xi_t < shape_lower - 1e-10)) {
    value <- 1e9
  } else {
    z_u    <- 1 + xi_t * (threshold - mu_t) / sig_t
    gumbel <- abs(xi_t) < xi_tol
    if (any(!is.finite(z_u)) || any(!gumbel & xi_t > 0 & z_u <= 0)) {
      value <- 1e9
    } else {
      tau    <- numeric(n)
      pos    <- !gumbel & z_u > 0
      tau[gumbel] <- exp(pmax(-500, -(threshold - mu_t[gumbel]) / sig_t[gumbel]))
      tau[pos]    <- z_u[pos]^(-1 / xi_t[pos])
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

# Cached objective for optim(). fn computes only the value, since optim()
# calls it at every line-search trial (many of which fall outside the
# support); gr computes value and gradient and caches both, so fn at an
# accepted point is free.
.objective_factory <- function(dm, y, threshold, pv, eps_l1 = 1e-10,
                               obs_per_year = 365.25, shape_lower = -Inf) {
  v_par <- NULL; v_val <- NULL
  g_par <- NULL; g_val <- NULL
  fn <- function(par) {
    if (!is.null(g_par) && identical(par, g_par)) return(g_val$value)
    if (!is.null(v_par) && identical(par, v_par)) return(v_val)
    core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = FALSE,
                       shape_lower = shape_lower)
    val <- if (core$value >= 1e9) 1e9 else
      core$value + .penalty_value_grad(par, pv, eps_l1)$value
    v_par <<- par; v_val <<- val
    val
  }
  gr <- function(par) {
    if (!is.null(g_par) && identical(par, g_par)) return(g_val$grad)
    core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = TRUE,
                       shape_lower = shape_lower)
    pen  <- .penalty_value_grad(par, pv, eps_l1)
    value <- if (core$value >= 1e9) 1e9 else core$value + pen$value
    g_val <<- list(value = value, grad = core$grad + pen$grad)
    g_par <<- par
    g_val$grad
  }
  list(fn = fn, gr = gr)
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

# Linear map for the lower bound on the shape parameter. With a shape design
# made of an intercept and at most one covariate x, xi_t >= shape_lower for
# every t holds iff it holds at the rows where x is smallest and largest, so
# v = A beta_xi (shape at those rows) turns the constraint into box bounds
# v >= shape_lower. Returns NULL when there is no bound or the design has
# more shape covariates (the bound is then enforced only as a barrier).
.shape_bound_map <- function(dm, shape_lower) {
  if (!is.finite(shape_lower)) return(NULL)
  X    <- dm$X_xi
  p_xi <- ncol(X)
  if (p_xi > 2L || any(X[, 1L] != 1)) return(NULL)
  rows <- if (p_xi == 1L) 1L else c(which.min(X[, 2L]), which.max(X[, 2L]))
  A <- X[rows, , drop = FALSE]
  if (abs(det(A)) < 1e-12) return(NULL)
  list(A = A, Ainv = solve(A),
       idx = ncol(dm$X_mu) + ncol(dm$X_sigma) + seq_len(p_xi),
       lower = shape_lower)
}

# KKT multipliers of the shape block under the bound. The gradient is mapped
# to the bounded coordinates v (h = Ainv' g); a coordinate at its bound must
# have h >= 0 (its multiplier), an interior one h = 0. A penalized shape
# covariate at zero contributes the subgradient lambda * alpha * s,
# s in [-1, 1]. The violation max_j viol_j(s) is convex and piecewise linear
# in s, so its minimum is attained at an endpoint, at a root of some h_j(s),
# or where two pieces cross; all such candidates are evaluated exactly (a
# grid fails when lambda * w is large, as with adaptive weights).
# Returns the violation, the multipliers h at the best s, and the bound flags.
.shape_kkt <- function(par, grad_total, pv, active_tol, bmap,
                       eps_l1 = 1e-10) {
  S <- bmap$idx
  v <- drop(bmap$A %*% par[S])
  at_bound <- v - bmap$lower < 1e-7
  l1   <- pv$lam[S] * pv$alp[S] * par[S] / sqrt(par[S]^2 + eps_l1)
  kink <- pv$lam[S] > 0 & pv$alp[S] > 0 & abs(par[S]) < active_tol
  base <- grad_total[S] - l1
  base[!kink] <- base[!kink] + l1[!kink]          # smooth part off the kink
  h0 <- drop(crossprod(bmap$Ainv, base))
  e  <- numeric(length(S))
  e[kink] <- pv$lam[S][kink] * pv$alp[S][kink]
  hc <- drop(crossprod(bmap$Ainv, e))
  cand <- 0
  if (any(kink)) {
    cand <- c(-1, 0, 1)
    nz <- abs(hc) > 0
    cand <- c(cand, -h0[nz] / hc[nz])
    m <- length(h0)
    for (i in seq_len(m)) for (j in seq_len(m)) if (i < j)
      for (sg in c(1, -1)) {
        den <- hc[i] - sg * hc[j]
        if (abs(den) > 0) cand <- c(cand, (sg * h0[j] - h0[i]) / den)
      }
    cand <- unique(pmin(1, pmax(-1, cand)))
  }
  viol <- function(s) {
    h <- h0 + s * hc
    max(ifelse(at_bound, pmax(0, -h), abs(h)))
  }
  vals <- vapply(cand, viol, numeric(1L))
  s <- cand[which.min(vals)]
  list(violation = min(vals), h = h0 + s * hc, at_bound = at_bound)
}

# KKT violation under the shape bound: .kkt_violation() outside the shape
# block, .shape_kkt() inside it.
.kkt_violation_bounded <- function(par, grad_total, pv, active_tol, bmap,
                                   eps_l1 = 1e-10) {
  if (is.null(bmap)) return(.kkt_violation(par, grad_total, pv, active_tol, eps_l1))
  if (any(!is.finite(grad_total))) return(Inf)
  S <- bmap$idx
  v <- drop(bmap$A %*% par[S])
  if (!any(v - bmap$lower < 1e-7))
    return(.kkt_violation(par, grad_total, pv, active_tol, eps_l1))
  other <- setdiff(seq_along(par), S)
  v_other <- if (length(other))
    .kkt_violation(par[other], grad_total[other],
                   list(lam = pv$lam[other], alp = pv$alp[other]),
                   active_tol, eps_l1) else 0
  max(v_other, .shape_kkt(par, grad_total, pv, active_tol, bmap, eps_l1)$violation)
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
#' @param shape_lower Lower bound for the shape parameter at every period;
#'   parameter vectors with some \eqn{\xi_t} below it return \code{1e9}.
#'   Default \code{-Inf} (no bound).
#'
#' @details The integrated intensity is
#'   \eqn{\Lambda = n_y^{-1}\sum_t [1 + \xi_t (u - \mu_t)/\sigma_t]_+^{-1/\xi_t}},
#'   so that, when \eqn{\xi_t < 0} and the upper endpoint
#'   \eqn{\mu_t - \sigma_t/\xi_t} lies below \eqn{u}, period \eqn{t}
#'   contributes zero (the intensity vanishes above the threshold).
#'
#' @return Scalar. The penalized negative log-likelihood value. Parameter
#'   vectors outside the point-process support (an exceedance beyond the
#'   upper endpoint, or a threshold below the lower endpoint) return
#'   \code{1e9}.
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
                    penalty_weights = NULL,
                    shape_lower = -Inf) {
  core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = FALSE,
                     shape_lower = shape_lower)
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
                    penalty_weights = NULL,
                    shape_lower = -Inf) {
  core <- .nhpp_eval(par, dm, y, threshold, obs_per_year, grad = TRUE,
                     shape_lower = shape_lower)
  pv <- .penalty_vectors(dm, lambda, alpha, pen_xi, penalty_weights)
  core$grad + .penalty_value_grad(par, pv, eps_l1)$grad
}
