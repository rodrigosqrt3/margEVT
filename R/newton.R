# =============================================================================
# newton.R
# Proximal Newton solver for the penalized point-process likelihood
# (fit_nhpp(solver = "newton")).
#
# At each iteration the negative log-likelihood is replaced by its
# second-order expansion (analytical Hessian, nhpp_hess_cpp()), the
# elastic-net penalty is kept exact, and the resulting penalized quadratic
# is minimized by cyclic coordinate descent with soft thresholding, as in
# glmnet. A backtracking line search on the exact objective guarantees
# descent. Penalized coefficients reach exact zeros, which the smoothed-L1
# BFGS solver only approaches slowly. Convergence is declared with the same
# KKT check as the BFGS solver.
#
#   - .nhpp_hessian() : value, gradient and Hessian in the coefficients
#   - .cd_quadratic() : coordinate descent for the penalized quadratic model
#   - .prox_newton()  : the solver
# =============================================================================

# Value, gradient and Hessian of the unpenalized negative log-likelihood in
# the coefficients, or NULL when the compiled kernel cannot be used (Gumbel
# limit, outside the support).
.nhpp_hessian <- function(par, dm, y, threshold, obs_per_year, shape_lower) {
  r <- nhpp_hess_cpp(as.numeric(par), dm$X_mu, dm$X_sigma, dm$X_xi,
                     as.numeric(y), as.numeric(threshold),
                     as.numeric(obs_per_year), as.numeric(shape_lower))
  if (r$status != 0L) return(NULL)
  Xm <- dm$X_mu; Xs <- dm$X_sigma; Xx <- dm$X_xi
  cw <- function(A, w, B) crossprod(A * w, B)
  Hmm <- cw(Xm, r$hmm, Xm); Hms <- cw(Xm, r$hms, Xs); Hmx <- cw(Xm, r$hmx, Xx)
  Hss <- cw(Xs, r$hss, Xs); Hsx <- cw(Xs, r$hsx, Xx); Hxx <- cw(Xx, r$hxx, Xx)
  H <- rbind(cbind(Hmm, Hms, Hmx),
             cbind(t(Hms), Hss, Hsx),
             cbind(t(Hmx), t(Hsx), Hxx))
  dimnames(H) <- NULL
  list(value = r$value, grad = r$grad, hessian = (H + t(H)) / 2)
}

# Exact elastic-net penalty sum_k lam_k (alp_k |b_k| + (1 - alp_k) b_k^2 / 2).
.penalty_exact <- function(b, pv)
  sum(pv$lam * (pv$alp * abs(b) + 0.5 * (1 - pv$alp) * b^2))

# Minimizes g'(b - theta) + (b - theta)' H (b - theta) / 2 + penalty(b)
# subject to b >= lower, by cyclic coordinate descent with soft thresholding.
.cd_quadratic <- function(theta, g, H, pv, lower, max_sweeps = 1000L,
                          tol = 1e-12) {
  p  <- length(theta)
  b  <- theta
  d  <- numeric(p)
  Hd <- numeric(p)
  l1 <- pv$lam * pv$alp
  l2 <- pv$lam * (1 - pv$alp)
  hd <- diag(H)
  for (sweep in seq_len(max_sweeps)) {
    moved <- 0
    for (j in seq_len(p)) {
      zj <- hd[j] * theta[j] - (g[j] + Hd[j] - hd[j] * d[j])
      bj <- if (l1[j] > 0) sign(zj) * max(0, abs(zj) - l1[j]) else zj
      bj <- max(bj / (hd[j] + l2[j]), lower[j])
      delta <- bj - b[j]
      if (delta != 0) {
        b[j]  <- bj
        d[j]  <- d[j] + delta
        Hd    <- Hd + H[, j] * delta
        moved <- max(moved, abs(delta))
      }
    }
    if (moved < tol * max(1, max(abs(b)))) break
  }
  b
}

# Proximal Newton for one penalty level. Returns a list like a run of
# .fit_at_lambda()'s optimizers (par, value, kkt), or NULL when the solver
# does not apply or fails, so that the caller falls back to BFGS.
# The shape bound is handled as a box when the shape design is an intercept
# only; with shape covariates under a finite bound the solver is not used.
.prox_newton <- function(init, dm, y, threshold, pv, obs_per_year,
                         shape_lower, bmap, tol_abs, active_tol,
                         max_iter = 200L) {
  p_mu <- ncol(dm$X_mu); p_sig <- ncol(dm$X_sigma); p_xi <- ncol(dm$X_xi)
  lower <- rep(-Inf, p_mu + p_sig + p_xi)
  if (is.finite(shape_lower)) {
    if (is.null(bmap) || p_xi != 1L) return(NULL)
    lower[bmap$idx] <- shape_lower
  }
  objective <- function(b) {
    v <- .nhpp_eval(b, dm, y, threshold, obs_per_year, grad = FALSE,
                    shape_lower = shape_lower)$value
    if (!is.finite(v) || v >= 1e9) Inf else v + .penalty_exact(b, pv)
  }
  par <- pmax(as.numeric(init), lower)
  f   <- objective(par)
  if (!is.finite(f)) return(NULL)
  for (it in seq_len(max_iter)) {
    h <- .nhpp_hessian(par, dm, y, threshold, obs_per_year, shape_lower)
    if (is.null(h)) return(NULL)
    g_tot <- h$grad + .penalty_value_grad(par, pv)$grad
    kkt <- .kkt_violation_bounded(par, g_tot, pv, active_tol, bmap)
    if (kkt <= tol_abs)
      return(list(par = stats::setNames(par, names(init)), value = f,
                  kkt = kkt, convergence = 0L, iterations = it - 1L))
    H  <- h$hessian
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    floor_ev <- 1e-8 * max(abs(ev), 1)
    if (min(ev) < floor_ev) diag(H) <- diag(H) + floor_ev - min(ev)
    b <- .cd_quadratic(par, h$grad, H, pv, lower)
    d <- b - par
    if (max(abs(d)) < 1e-14) return(NULL)
    decrease <- sum(h$grad * d) + .penalty_exact(b, pv) -
      .penalty_exact(par, pv)
    step <- 1
    repeat {
      cand <- par + step * d
      f_c  <- objective(cand)
      if (is.finite(f_c) && f_c <= f + 1e-4 * step * min(decrease, 0)) break
      step <- step / 2
      if (step < 1e-12) return(NULL)
    }
    par <- cand
    f   <- f_c
  }
  NULL
}
