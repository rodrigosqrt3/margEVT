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
