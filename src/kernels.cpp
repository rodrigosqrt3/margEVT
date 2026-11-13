// Compiled kernels for the two hot loops of margEVT:
//
//   nhpp_eval_cpp()  negative log-likelihood of the point-process model and
//                    its gradient, in one pass over the observations
//                    (R wrapper: .nhpp_eval() in likelihood.R);
//   marg_paths_cpp() pathwise annual non-exceedance probabilities and their
//                    derivatives in the level (R wrapper: .marg_paths() in
//                    marginalize.R).
//
// Both reproduce the R definitions they replace (.nhpp_eval_ref(),
// .tail_measure_at_level(), .tail_measure_derivative()) to rounding error.

#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Linear predictor X %*% b for a column-major matrix, skipping zero
// coefficients (most of them after selection).
static void lin_pred(const NumericMatrix &X, const double *b,
                     std::vector<double> &out) {
  const R_xlen_t n = X.nrow();
  const int p = X.ncol();
  out.assign(n, 0.0);
  const double *x = X.begin();
  for (int j = 0; j < p; ++j) {
    const double bj = b[j];
    if (bj == 0.0) continue;
    const double *col = x + (R_xlen_t)j * n;
    for (R_xlen_t i = 0; i < n; ++i) out[i] += col[i] * bj;
  }
}

static void cross_prod(const NumericMatrix &X, const std::vector<double> &g,
                       double *out) {
  const R_xlen_t n = X.nrow();
  const int p = X.ncol();
  const double *x = X.begin();
  for (int j = 0; j < p; ++j) {
    const double *col = x + (R_xlen_t)j * n;
    double acc = 0.0;
    for (R_xlen_t i = 0; i < n; ++i) acc += col[i] * g[i];
    out[j] = acc;
  }
}

// status: 0 = value (and gradient) computed; 1 = outside the model support
// (value 1e9); 2 = defer to the R reference implementation (a non-finite or
// Gumbel-limit shape, |xi| < 1e-6, at some period).
// [[Rcpp::export]]
List nhpp_eval_cpp(NumericVector par, NumericMatrix Xm, NumericMatrix Xs,
                   NumericMatrix Xx, NumericVector y, double u, double ny,
                   bool grad, double shape_lower) {
  const R_xlen_t n = y.size();
  const int pm = Xm.ncol(), ps = Xs.ncol(), px = Xx.ncol();
  const double *b = par.begin();
  std::vector<double> mu, ls, xi;
  lin_pred(Xm, b, mu);
  lin_pred(Xs, b + pm, ls);
  lin_pred(Xx, b + pm + ps, xi);

  for (R_xlen_t i = 0; i < n; ++i)
    if (!std::isfinite(xi[i]) || std::fabs(xi[i]) < 1e-6)
      return List::create(_["status"] = 2);

  std::vector<double> gm, gs, gx;
  if (grad) { gm.assign(n, 0.0); gs.assign(n, 0.0); gx.assign(n, 0.0); }
  const List outside = List::create(_["status"] = 1, _["value"] = 1e9);

  double lam = 0.0, logf = 0.0;
  const double *yy = y.begin();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double s = std::exp(ls[i]), k = xi[i];
    if (!std::isfinite(mu[i]) || !std::isfinite(s) || s <= 0.0 ||
        k < shape_lower - 1e-10) return outside;
    const double wu = (u - mu[i]) / s, zu = 1.0 + k * wu;
    if (!std::isfinite(zu) || (k > 0.0 && zu <= 0.0)) return outside;
    if (zu > 0.0) {
      const double lz = std::log(zu), tau = std::exp(-lz / k);
      lam += tau;
      if (grad) {
        const double tz = tau / zu;
        gm[i] = tz / s / ny;
        gs[i] = tz * wu / ny;
        gx[i] = tau * (lz / (k * k) - wu / (zu * k)) / ny;
      }
    }
    if (yy[i] > u) {
      const double we = (yy[i] - mu[i]) / s, ze = 1.0 + k * we;
      if (!std::isfinite(ze) || ze <= 0.0) return outside;
      const double lze = std::log(ze);
      logf += -ls[i] - (1.0 + 1.0 / k) * lze;
      if (grad) {
        const double a = (1.0 + k) / ze;
        gm[i] -= a / s;
        gs[i] += 1.0 - a * we;
        gx[i] += -lze / (k * k) + (1.0 + 1.0 / k) * we / ze;
      }
    }
  }
  const double value = lam / ny - logf;
  if (!grad) return List::create(_["status"] = 0, _["value"] = value);

  NumericVector g(pm + ps + px);
  cross_prod(Xm, gm, g.begin());
  cross_prod(Xs, gs, g.begin() + pm);
  cross_prod(Xx, gx, g.begin() + pm + ps);
  return List::create(_["status"] = 0, _["value"] = value, _["grad"] = g);
}

// Value, gradient and the six per-period weights of the Hessian of the
// negative log-likelihood with respect to the linear predictors
// (mu, log sigma, xi): h_mm, h_ms, h_mx, h_ss, h_sx, h_xx. The Hessian in the
// coefficients is then sum_t J_t' h_t J_t (assembled in R). Only for the
// interior case: status 2 (Gumbel limit or non-finite shape) and 1 (outside
// the support) are returned as in nhpp_eval_cpp().
//
// With w = (u - mu)/sigma, z = 1 + xi w, L = log tau = -log(z)/xi, the
// compensator of period t contributes tau (L_a L_b + L_ab) / n_y, and an
// exceedance with v = (y - mu)/sigma, s = 1 + xi v, M = log s and
// c = 1 + 1/xi contributes the second derivatives of log sigma + c M.
// [[Rcpp::export]]
List nhpp_hess_cpp(NumericVector par, NumericMatrix Xm, NumericMatrix Xs,
                   NumericMatrix Xx, NumericVector y, double u, double ny,
                   double shape_lower) {
  const R_xlen_t n = y.size();
  const int pm = Xm.ncol(), ps = Xs.ncol(), px = Xx.ncol();
  const double *b = par.begin();
  std::vector<double> mu, ls, xi;
  lin_pred(Xm, b, mu);
  lin_pred(Xs, b + pm, ls);
  lin_pred(Xx, b + pm + ps, xi);
  for (R_xlen_t i = 0; i < n; ++i)
    if (!std::isfinite(xi[i]) || std::fabs(xi[i]) < 1e-6)
      return List::create(_["status"] = 2);

  std::vector<double> gm(n, 0.0), gs(n, 0.0), gx(n, 0.0);
  NumericVector hmm(n), hms(n), hmx(n), hss(n), hsx(n), hxx(n);
  const List outside = List::create(_["status"] = 1, _["value"] = 1e9);
  double lam = 0.0, logf = 0.0;
  const double *yy = y.begin();
  for (R_xlen_t i = 0; i < n; ++i) {
    const double s = std::exp(ls[i]), k = xi[i];
    if (!std::isfinite(mu[i]) || !std::isfinite(s) || s <= 0.0 ||
        k < shape_lower - 1e-10) return outside;
    const double w = (u - mu[i]) / s, z = 1.0 + k * w;
    if (!std::isfinite(z) || (k > 0.0 && z <= 0.0)) return outside;
    if (z > 0.0) {
      const double lz = std::log(z), tau = std::exp(-lz / k), t = tau / ny;
      lam += tau;
      const double Lm = 1.0 / (s * z), Ls = w / z,
        Lx = lz / (k * k) - w / (k * z);
      const double Lmm = k / (s * s * z * z),
        Lms = -1.0 / (s * z) + k * w / (s * z * z),
        Lmx = -w / (s * z * z),
        Lss = -w / z + k * w * w / (z * z),
        Lsx = -w * w / (z * z),
        Lxx = -2.0 * lz / (k * k * k) + 2.0 * w / (k * k * z) +
          w * w / (k * z * z);
      gm[i] = t * Lm; gs[i] = t * Ls; gx[i] = t * Lx;
      hmm[i] = t * (Lm * Lm + Lmm); hms[i] = t * (Lm * Ls + Lms);
      hmx[i] = t * (Lm * Lx + Lmx); hss[i] = t * (Ls * Ls + Lss);
      hsx[i] = t * (Ls * Lx + Lsx); hxx[i] = t * (Lx * Lx + Lxx);
    }
    if (yy[i] > u) {
      const double v = (yy[i] - mu[i]) / s, se = 1.0 + k * v;
      if (!std::isfinite(se) || se <= 0.0) return outside;
      const double M = std::log(se), c = 1.0 + 1.0 / k, k2 = k * k;
      logf += -ls[i] - c * M;
      const double Mm = -k / (s * se), Ms = -k * v / se, Mx = v / se;
      const double Mmm = -k2 / (s * s * se * se), Mms = k / (s * se * se),
        Mmx = -1.0 / (s * se * se), Mss = k * v / (se * se),
        Msx = -v / (se * se), Mxx = -v * v / (se * se);
      gm[i] += c * Mm; gs[i] += 1.0 + c * Ms; gx[i] += -M / k2 + c * Mx;
      hmm[i] += c * Mmm; hms[i] += c * Mms; hss[i] += c * Mss;
      hmx[i] += c * Mmx - Mm / k2;
      hsx[i] += c * Msx - Ms / k2;
      hxx[i] += c * Mxx - 2.0 * Mx / k2 + 2.0 * M / (k2 * k);
    }
  }
  NumericVector g(pm + ps + px);
  cross_prod(Xm, gm, g.begin());
  cross_prod(Xs, gs, g.begin() + pm);
  cross_prod(Xx, gx, g.begin() + pm + ps);
  return List::create(_["status"] = 0, _["value"] = lam / ny - logf,
                      _["grad"] = g, _["hmm"] = hmm, _["hms"] = hms,
                      _["hmx"] = hmx, _["hss"] = hss, _["hsx"] = hsx,
                      _["hxx"] = hxx);
}

// For a stack of paths of n_rows consecutive periods, the probability that
// level z is not exceeded along each path, G_b = exp(-sum_t tau_t(z) / n_obs),
// and (deriv = TRUE) its derivative in z, dG_b = -G_b sum_t tau_t'(z) / n_obs
// (zero when G_b = 0). tau_t(z) = [1 + xi_t (z - mu_t) / sigma_t]^(-1/xi_t):
// infinite below the lower endpoint (xi_t > 0), zero above the upper
// endpoint (xi_t < 0), exp(-(z - mu_t) / sigma_t) in the Gumbel limit.
// Parameters are validated in R (finite, positive scales).
// [[Rcpp::export]]
List marg_paths_cpp(double z, NumericVector mu, NumericVector sigma,
                    NumericVector xi, int n_rows, double n_obs, bool deriv) {
  const R_xlen_t n = mu.size();
  const R_xlen_t B = n / n_rows;
  NumericVector G(B), dG(deriv ? B : 0);
  const double *m = mu.begin(), *s = sigma.begin(), *k = xi.begin();
  for (R_xlen_t b = 0; b < B; ++b) {
    double S = 0.0, D = 0.0;
    bool infinite = false;
    const R_xlen_t lo = b * n_rows, hi = lo + n_rows;
    for (R_xlen_t i = lo; i < hi; ++i) {
      const double w = (z - m[i]) / s[i];
      if (std::fabs(k[i]) < 1e-6) {
        const double t = std::exp(std::max(-500.0, -w));
        S += t;
        D -= t / s[i];
        continue;
      }
      const double inner = 1.0 + k[i] * w;
      if (inner > 0.0) {
        const double li = std::log(inner);
        const double t = std::exp(-li / k[i]);
        S += t;
        D -= t / inner / s[i];
      } else if (k[i] > 0.0) {
        infinite = true;
      }
    }
    if (infinite || !std::isfinite(S)) {
      G[b] = 0.0;
      if (deriv) dG[b] = 0.0;
    } else {
      G[b] = std::exp(-S / n_obs);
      if (deriv) dG[b] = G[b] > 0.0 ? -G[b] * D / n_obs : 0.0;
    }
  }
  return List::create(_["G"] = G, _["dG"] = dG);
}
