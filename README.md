# margEVT: Regularized Point Processes and Stochastic Marginalization for Extremes

[![CRAN status](https://www.r-pkg.org/badges/version/margEVT)](https://CRAN.R-project.org/package=margEVT) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/margEVT/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/margEVT/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/margEVT/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/margEVT)

`margEVT` fits non-stationary extreme value models in which climate indices
and other covariates drive the parameters, selects those covariates by
penalized likelihood, and turns the fitted model into return levels and
risk measures that integrate the covariates out, instead of fixing them at
an arbitrary value.

## What it does

| Task | Functions |
|---|---|
| Fit a non-homogeneous Poisson process (NHPP) for threshold exceedances, with LASSO, adaptive LASSO or elastic-net selection of covariates, BIC along the regularization path and an optional unpenalized refit | `fit_nhpp()`, `refit_nhpp()`, `seasonal_harmonics()` |
| Inspect and diagnose the fit | `summary()`, `coef()`, `predict()`, `fitted()`, `residuals()`, `plot()`, `logLik()`, `AIC()`, `BIC()`, `vcov()`, `confint()`, `simulate()`, `update()` |
| Handle clustered exceedances | `extremal_index()`, `decluster()` |
| Return levels: conditional on covariate scenarios (A) or marginal over a VAR generator (B) or resampled years (C) | `marginalize()`, `fit_var_generator()`, `simulate_covariates()`, `rl_table()`, `annual_prob()` |
| Risk over a planning horizon of several years, stationary or not | `horizon_risk()`, `design_life_level()` |
| Uncertainty and validation | `bootstrap_nhpp()`, `bootstrap_rl()`, `bootstrap_coef()`, `backtest()` |

## Installation

From CRAN:

```r
install.packages("margEVT")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("rodrigosqrt3/margEVT")
```

## Quick start

```r
library(margEVT)

# 40 years of weekly data: annual cycle and an autocorrelated index x
set.seed(2024)
n_y <- 52; n <- 40 * n_y
df <- cbind(
  data.frame(year = rep(1981:2020, each = n_y),
             x = as.numeric(arima.sim(list(ar = 0.8), n)) * 0.6,
             z = rnorm(n)),                       # z is irrelevant
  seasonal_harmonics(rep(seq_len(n_y), 40), period = n_y)
)
df$y <- 20 + 3 * df$x + 4 * df$cos1 + 3 * rexp(n)
u <- quantile(df$y, 0.9)
df$y <- decluster(df$y, u)                        # one value per cluster

# Adaptive LASSO with BIC selection; the harmonics are not penalized
fit <- fit_nhpp(df, threshold = u,
                loc_vars = c("cos1", "sen1", "x", "z"),
                scale_vars = c("x", "z"), free_vars = c("cos1", "sen1"),
                penalty = "lasso", lambda = "bic",
                penalty_weights = "adaptive", refit = TRUE,
                obs_per_year = n_y, verbose = FALSE)
summary(fit)
plot(fit, type = "qq")

# Return levels: at fixed values of x (A), over a VAR generator (B),
# over resampled observed years (C)
gen  <- fit_var_generator(fit, df, vars = "x")
sims <- simulate_covariates(gen, n_mc = 300, seed = 1)
rl <- marginalize(fit, df, TRs = c(10, 50, 100),
                  approaches = c("A", "B", "C"),
                  scenarios = list(neutral = list(x = 0)),
                  mc_sample = sims, n_boot = 300, seed = 1)
rl_table(rl)

# Probability of exceeding 45 at least once in the next 10 years
traj <- replicate(200, lapply(sample(1981:2020, 10, replace = TRUE),
                              function(yr) df[df$year == yr, "x", drop = FALSE]),
                  simplify = FALSE)
horizon_risk(fit, z = 45, trajectories = traj)
```

The vignette, `vignette("margEVT")`, goes through the same analysis step by
step, with diagnostics, bootstrap intervals and walk-forward validation.

## The model in brief

Exceedances of a threshold $u$ follow a point process whose intensity is
parameterized by the generalized extreme value (GEV) parameters of the
annual maximum,

$$\Lambda_t\{(y, \infty)\} = \frac{1}{n_y}\left[1 + \xi_t\,\frac{y - \mu_t}{\sigma_t}\right]_+^{-1/\xi_t}, \qquad y > u,$$

with $\mu_t$, $\log\sigma_t$ and $\xi_t$ linear in the covariates and
$n_y$ observations per year. The coefficients minimize the negative
log-likelihood plus a weighted elastic-net penalty; intercepts and seasonal
harmonics are left unpenalized, the shape is kept at $\xi_t \ge -1/2$
(the region where maximum likelihood is regular; Smith, 1985), and every
fit is checked against the Karush-Kuhn-Tucker conditions.

Given one year of covariates $\mathbf{v}$, the annual non-exceedance
probability is $G(z \mid \mathbf{v}) = \exp\{-\sum_t \Lambda_t(z)\}$. A
level computed at a fixed $\mathbf{v}$ is conditional on that state. The
marginal $T$-year return level instead averages over a law $\Pi$ for whole
years of covariates,

$$G_\Pi(z_T) = \mathrm{E}_\Pi\{G(z_T \mid \mathbf{V})\} = 1 - 1/T,$$

with $\Pi$ given by a stable VAR(p) fitted to the deseasonalized covariates
(approach B) or by the observed years (approach C). Over a horizon of $L$
years with trajectories $\mathbf{v}_{1s}, \ldots, \mathbf{v}_{Ls}$, the
probability of no exceedance is $S^{-1}\sum_s \prod_l G(z \mid
\mathbf{v}_{ls})$, which gives the horizon risk and the design life level
(Rootzén and Katz, 2013).

## Citation

Run `citation("margEVT")` in R, or cite:

> Villa, R. F. (2026). *Regularized Poisson Point Process Models and
> Stochastic Marginalization for Return-Level Inference under
> Covariate-Driven Extremes*. Master's dissertation, Instituto de Matemática
> e Estatística, Universidade Federal do Rio Grande do Sul, Porto Alegre,
> Brazil. Advisor: Flávio A. Ziegelmann.

## References

Rootzén, H. and Katz, R. W. (2013). Design life level: quantifying risk in a
changing climate. *Water Resources Research*, 49(9), 5964–5972.

Smith, R. L. (1985). Maximum likelihood estimation in a class of nonregular
cases. *Biometrika*, 72(1), 67–90.
