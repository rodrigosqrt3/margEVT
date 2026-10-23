# margEVT 0.5.0

## Breaking changes

- The point-process likelihood now uses the `[.]_+` convention of the
  integrated intensity: when `xi_t < 0` and the upper
  endpoint `mu_t - sigma_t / xi_t` lies below the threshold, period `t`
  contributes zero to `Lambda`. Previously such parameter vectors were
  rejected (`1e9`) while the gradient already treated them as zero, so fits
  with a negative shape could stop at an artificial boundary without being
  stationary. Values remain `1e9` only outside the model support (an
  exceedance beyond the upper endpoint, or a threshold below the lower
  endpoint when `xi_t > 0`).
- New `fit_nhpp(shape_lower = -0.5)`: the shape parameter is restricted to
  `xi_t >= shape_lower` (Smith, 1985, regular region; below `-1` the
  likelihood is unbounded). For a shape design with an intercept and at
  most one covariate the bound is handled exactly, with KKT conditions and
  multipliers; otherwise it acts as a barrier. Use `shape_lower = -Inf` for
  the previous unbounded behaviour. Fitted objects store `shape_lower` and
  `shape_at_bound`, and `print()` reports an active bound.
- `pp_nllh()` and `pp_grad()` gain `shape_lower` (default `-Inf`).
- Bootstrap and backtest refits inherit `shape_lower` from the fit.
- `fit_nhpp()` gains `scaling`, with default `"auto"` (see Performance).
  Estimates change within the optimizer tolerance, and selected models were
  unchanged in all tests; `scaling = "none"` reproduces earlier versions.

## Estimation

- New active-set solver for the shape bound: bounded shape coordinates at
  the bound are fixed there and the remaining coordinates are optimized by
  BFGS, with the active set updated from the multiplier signs. L-BFGS-B,
  used before, stalled near the smoothed L1 kink: on a daily temperature
  refit at the bound it ran ~250 s without passing the KKT check, against
  ~20 s and a certified (lower) optimum now. Fits that start at the bound
  (along a regularization path, or in bootstrap and backtest refits once the
  bound is active) go straight to this solver; L-BFGS-B remains a fallback.
- KKT check under the shape bound: the subgradient of a penalized shape
  covariate at zero is now chosen exactly (the violation is piecewise linear
  in it) instead of on a grid. With adaptive weights (lambda * w of order
  1e5) the grid missed the admissible interval and optimal fits were
  reported as non-converged, leaving gaps in the regularization path (e.g.
  27 of 56 path points for a daily temperature fit).
- When the KKT check fails, BFGS is first restarted from its solution (fresh
  Hessian approximation) before falling back to L-BFGS-B; this polishes fits
  that stop just short of the tolerance and removes most failed path points.

## Performance

- The likelihood is evaluated in a single pass (`z_u`, `log z_u` and `tau`
  once; `z_u^(-1/xi - 1)` as `tau / z_u`), about 4 times faster per call.
- The objective passed to `optim()` computes the gradient only when it is
  requested: line-search trials, many of which fall outside the support,
  now cost one value evaluation. Previously every call also computed the
  gradient, and support violations went through the slow path.
- Together: an unpenalized fit on 1000 observations went from ~60 ms to
  ~19 ms, and an adaptive-lasso fit with the BIC path from ~5.0 s to ~1.4 s,
  with the same estimates (differences of order 1e-8). The previous
  implementation is kept as `.nhpp_eval_ref()` for support violations with
  gradient and for Gumbel limits (`|xi| < 1e-6`).
- `fit_nhpp()` sets `options(matprod = "blas")` while it runs (restored on
  exit). Its covariates are checked to be finite, so the NaN scan that the
  default runs over the whole design matrix before each product is
  redundant; results are bitwise identical. On a daily-size adaptive fit
  (16436 days, 52 parameters) this cut the time from 182 s to 146 s.
- New argument `scaling` in `fit_nhpp()`. With `"auto"` (default) the
  optimizer receives the parameter scales `1 / sqrt(|H_kk|)`, from the
  diagonal curvature of the unpenalized likelihood at the starting values
  (computed once per call, reused along the regularization path and stored
  in the fit, so bootstrap and backtest refits use the same setting).
  Unscaled BFGS evaluated the objective about 9 times per gradient, cutting
  back steps of the wrong size. On the daily-size fit above: 147 s -> 77 s,
  same lambda and active set, lower objective (negative log-likelihood
  1093.94 vs 1094.09) and a KKT violation 28 times smaller (6.8e-6 vs
  1.9e-4), with all 54 path points certified (53 before). Estimates change
  within the optimizer tolerance; along flat directions of the likelihood
  the change can reach a few percent of a coefficient (2.6% for one
  location slope in that example), always towards a lower objective.
  `scaling = "none"` reproduces earlier versions bitwise.
- Benchmark of the JSS article (unpenalized fit, n = 1000): 62 ms -> 17 ms
  (extRemes 47 ms, ismev 176 ms on the same machine).

## Bootstrap

- `bootstrap_rl()`, `bootstrap_coef()` and `bootstrap_nhpp()` gain
  `occurrence`: `"fixed"` (default, as before) keeps the observed number and
  timing of exceedances; `"simulated"` also draws them from the fitted point
  process (probability `1 - exp(-tau_t / n_y)` per period), so the
  uncertainty about the exceedance rate enters the intervals.
- New `bootstrap_nhpp()`: return-level and coefficient intervals from one set
  of refits. New `refit_nhpp()` and `annual_prob()` (exported helpers).
- Arguments are evaluated before replicates are sent to worker processes
  (lazily evaluated arguments failed on PSOCK workers).

## Backtesting

- `backtest()` gains `approaches`: `"C"` (default; resampled training years,
  as before) and/or `"B"`, a VAR generator fitted to the training years only
  (`generator_args` are passed to `fit_var_generator()`). Both approaches are
  evaluated on the same window refit; results carry `F_ann_<a>` and
  `exc_<a>_T<TR>` columns, `binom_tests` and `ks_test` an `approach` column.
  The legacy columns `F_ann` and `exc_T<TR>` repeat the first approach.
- `backtest()` gains `decluster`, a function applied to the training-year
  series in each window, so declustering no longer uses validation data.

## Minor

- `build_cov_annual()` no longer warns about a covariate missing from
  `cov_vals` when all its operational coefficients are zero (filling it with
  0 is then exact). This happens whenever the VAR generator, which models
  only the active covariates, feeds trajectories to a fit with inactive
  candidates.

# margEVT 0.4.0

## Breaking changes

- `fit_nhpp(lambda_scaling = )` now defaults to `"common"` (one scalar
  penalty for all parameter blocks), the convention used in the
  dissertation and articles. Use `lambda_scaling = "gradient"` for the
  previous score-calibrated block penalties.
- With `lambda = "bic"`, the returned model is now the regularization-path
  solution at the selected penalty, so the reported fit is exactly the one
  scored by the BIC.
- Convergence is now certified by the Karush-Kuhn-Tucker conditions of the
  penalized problem (`grad_tol`, scaled by the number of exceedances)
  rather than by the optimizer's return code. Fits that previously reported
  convergence may now warn, and vice versa.
- `summary()` returns an object of class `summary.nhpp_fit` (printed as
  before) instead of printing directly.
- The `vars` dependency was removed. VAR(p) generators are fitted by
  least squares with Schwarz order selection on a common sample, which
  reproduces `vars::VARselect()` and `vars::VAR()`; a single covariate is a
  VAR(p) with K = 1 rather than a separate AR model.
- Random-number streams changed: bootstrap replicates use independent
  L'Ecuyer-CMRG streams and trajectories are simulated jointly, so results
  for a given seed differ from 0.3.0.

## Estimation

- Adaptive LASSO: `penalty_weights = "adaptive"` computes weights
  `1 / |theta|^adaptive_gamma` from an unpenalized fit (ridge fallback);
  numeric coordinate weights are also accepted.
- New `fit_nhpp()` arguments `lambda_grid`, `init`, and `grad_tol`. The
  optimizer falls back from BFGS to L-BFGS-B and a restarted BFGS when the
  KKT check fails.
- Fitted objects store the penalty weights, the KKT violation, the
  evaluated regularization path (`fit$path`), and the observed series.
- The likelihood, gradient, and penalty are evaluated once per parameter
  vector and shared between objective and gradient calls.

## Methods and diagnostics

- New methods `logLik()`, `nobs()`, `predict()`, `residuals()`, and
  `coef()` for `nhpp_fit`; `AIC()` and `BIC()` work through `logLik()`.
- `residuals()` returns time-rescaling interarrival residuals and
  magnitude residuals, both standard exponential under the model.
- `plot.nhpp_fit()` gains `type = "path"` (BIC curve and coefficient
  paths) and `type = "qq"` (exponential Q-Q plots with exact pointwise
  bands).
- `n_exceedances()` uses the stored series when `y` is omitted.

## Marginalization and covariates

- New `seasonal_harmonics()`; harmonics of any order (`cos<k>`, `sin<k>`,
  `sen<k>`) are recognized and rebuilt automatically.
- `marginalize()` returns an `nhpp_rl` object with Monte Carlo standard
  errors of the return-level roots (`se_mc`), accepts `fixed` covariate
  values (e.g. `list(trend = 1)`), fills inactive design columns silently,
  and has a `plot()` method.
- `fit_var_generator()` gains `n_harmonics`, excludes deterministic
  variables from the VAR, and supports derived covariates (lags, moving
  averages) through `derive` and `derive_lead`. `simulate_covariates()`
  simulates all trajectories at once and gains `seed`.

## Bootstrap and walk-forward assessment

- Bootstrap refits start from the original estimate and keep the original
  penalty weights, gradient tolerance, and scaling. Results do not depend
  on the new `cores` argument; replicate matrices are returned as the
  attribute `"replicates"`.
- `backtest()` gains `trend_col` (trend renormalized within each training
  window), `fixed`, `reselect_lambda`, and `cores`. Adaptive weights are
  recomputed inside each training window, and covariates inactive in a
  window refit are filled without warnings.

## Documentation

- New vignette, examples for every exported function, and an updated
  citation file.

# margEVT 0.3.0

## Estimation consistency and reproducibility

- Added `lambda_scaling` to `fit_nhpp()`. Automatic BIC selection can retain
  the score-calibrated block penalties introduced in 0.2.0 or use a common
  scalar penalty to reproduce the estimation convention used in the
  dissertation analyses.
- Corrected automatic BIC selection so that both the unpenalized likelihood
  and the effective parameter count are evaluated at the thresholded
  operational estimator. Smooth optimizer solutions, rather than thresholded
  vectors, are retained as warm starts along the regularization path.
- Fitted objects now store both the smooth optimizer solution (`par`) and the
  operational activity-thresholded solution (`par_oper`). Downstream parameter
  prediction uses the operational estimator by default, while
  `coef(fit, operational = FALSE)` preserves access to the smooth solution.
- Coefficient bootstrap replicates now use the operational estimator, allowing
  percentile distributions to retain an atom at zero.
- Univariate covariate generators now fit a genuine stable AR model selected
  by BIC instead of padding the series with an artificial noise coordinate.
- Covariate simulation now supports exact Gaussian stationary-state
  initialization, while retaining zero-state burn-in as an explicit
  sensitivity option.

## Validation and defensive checks

- Return-level integration now rejects non-finite fitted parameters and
  non-positive scales instead of silently omitting invalid contributions.
- Parametric and empirical marginalization verify all direct and interaction
  input covariates before constructing annual trajectories.
- Prediction matrices are rebuilt directly from the stored design columns,
  avoiding formula reinterpretation of non-syntactic or interaction-like
  column names.
- Package metadata and the dissertation citation were synchronized for the
  0.3.0 development line.

# margEVT 0.2.0

## Consistency and temporal resolution

- Added `active_tol` to `fit_nhpp()` and stored it in fitted objects. The same
  tolerance is now used by automatic BIC selection, `print()`, `summary()`,
  `bic_nhpp()`, `active_covariates()`, bootstrap refits, and downstream
  covariate generators.
- `build_cov_annual()`, `fit_var_generator()`, `simulate_covariates()`,
  `marginalize()`, and `backtest()` now inherit the temporal resolution from
  the fitted model unless the user supplies an explicit override.
- Added frequency-aware behavior for weekly and other non-daily analyses.
- `fit_var_generator()` now verifies the VAR root condition, stores the root
  moduli and spectral radius, and refuses to simulate from an unstable fit.
- `simulate_covariates()` now generates independent annual Monte Carlo paths,
  each with its own burn-in, instead of splitting one dependent long path into
  nominal annual replicates.
- Corrected GEV endpoint handling in marginalized probabilities: levels below
  the lower endpoint when the shape is positive now have zero non-exceedance
  probability, while levels above a finite upper endpoint retain probability
  one conditionally on that trajectory.

## Inference and documentation

- Backtesting now reports `not_rejected` alongside the legacy `calibrated`
  field and uses language appropriate for nominal binomial comparisons and
  exploratory tail-conditional PIT diagnostics.
- Validation years with insufficient observed coverage are now excluded from
  binomial denominators instead of being silently counted as non-exceedances.
- Approach C now uses a locally preserved random seed, preventing its annual
  block sampling from resetting the outer bootstrap stream and repeating
  subsequent parametric bootstrap draws.
- `bootstrap_rl()` now reports the original fitted-model return level in
  `RL_est`; previous versions reported the mean of successful bootstrap roots.
- Clarified that the bootstrap is conditional on the observed occurrence
  pattern and selected model specification.
- Corrected the documentation of empirical annual-block marginalization and
  the event-based BIC convention.
- Added strict validation for block-specific penalty controls, return periods,
  Monte Carlo samples, and covariate trajectory lengths to prevent silent
  recycling or malformed fits.
- Return-level root finding no longer expands below the fitted threshold. It
  expands only the upper bracket and returns `NA` when the requested quantile
  lies outside the point-process tail domain.

# margEVT 0.1.1

- Added `plot.nhpp_fit()` for cumulative-intensity and fitted-parameter plots.
- Added package URL and bug-report metadata.
