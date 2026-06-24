# margEVT: Regularized Point Processes and Stochastic Marginalization for Return Level Inference

[![CRAN status](https://www.r-pkg.org/badges/version/margEVT)](https://CRAN.R-project.org/package=margEVT) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/margEVT/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/margEVT/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/margEVT/branch/main/graph/badge.svg)](https://codecov.io/gh/rodrigosqrt3/margEVT)

`margEVT` is an R package developed to conduct non-stationary extreme value analysis under covariate-driven regimes. The package implements the statistical framework developed in Villa (2026) under the supervision of Prof. Dr. Flavio Ziegelmann, coupling a covariate-driven Non-Homogeneous Poisson Process (NHPP) with an Elastic-Net penalized maximum likelihood estimation framework and stochastically marginalized return-level estimation.

---

## 1. Mathematical Framework

### 1.1 Non-Homogeneous Poisson Process (NHPP) Representation
Exceedances of an extreme threshold $u$ over an observational domain are modeled as a Non-Homogeneous Poisson Process on the product space $\mathcal{S} = [0, 1] \times (u, \infty)$. The time-varying intensity function $\lambda(t, y)$ is parameterized using the location $\mu(t)$, scale $\sigma(t) > 0$, and shape $\xi(t)$ parameters of the Generalized Extreme Value (GEV) distribution:

$$\lambda(t, y) = \frac{1}{\sigma(t)} \left[ 1 + \xi(t) \left( \frac{y - \mu(t)}{\sigma(t)} \right) \right]_{+}^{-1/\xi(t)-1}$$

where $[a]_+ = \max(a, 0)$. 

To capture non-stationarity, the parameters are modeled as linear combinations of covariate vectors $\mathbf{x}_t$:

$$\mu(t) = \mu_0 + \sum_{k=1}^{K} \mu_k x_{k,t}$$

$$\log(\sigma(t)) = \phi_0 + \sum_{j=1}^{J} \phi_j x_{j,t}$$

$$\xi(t) = \xi_0 + \sum_{l=1}^{L} \xi_l x_{l,t}$$

The logarithmic parameterization of $\sigma(t) = \exp(\phi(t))$ guarantees positivity.

### 1.2 Elastic-Net Regularization and Optimization
In high-dimensional covariate settings, parameter estimation is conducted via penalized maximum likelihood. The objective function minimizes the negative log-likelihood $\ell(\boldsymbol{\theta})$ subject to an Elastic-Net penalty:

$$Q(\boldsymbol{\theta}; \lambda, \alpha) = -\ell(\boldsymbol{\theta}) + \lambda \sum_{k \in \mathcal{P}} \left[ \alpha \sqrt{\theta_k^2 + \varepsilon} + \left( \frac{1 - \alpha}{2} \right) \theta_k^2 \right]$$

where:
- $\lambda \ge 0$ controls the overall regularization strength.
- $\alpha \in [0, 1]$ controls the balance between the $\ell_1$ (LASSO) and $\ell_2$ (Ridge) components (setting $\alpha = 1$ recovers the pure LASSO penalty).
- $\varepsilon > 0$ is a small smoothing parameter ensuring differentiability at the origin.
- $\mathcal{P}$ is the index set of penalized parameters. The baseline intercepts ($\mu_0, \phi_0, \xi_0$) and deterministic seasonal Fourier harmonics are kept unpenalized by construction ($\mathcal{P} \cap \text{unpenalized} = \emptyset$).

The optimization is solved via the quasi-Newton BFGS or L-BFGS-B algorithm supplied with exact analytical gradients:

$$\nabla_{\boldsymbol{\theta}} Q(\boldsymbol{\theta}; \lambda, \alpha) = -\nabla_{\boldsymbol{\theta}} \ell(\boldsymbol{\theta}) + \nabla_{\boldsymbol{\theta}} P_{\lambda,\alpha}(\boldsymbol{\theta})$$

The optimal regularization path parameter $\lambda^{\ast}$ is selected by minimizing the Bayesian Information Criterion (BIC):

$$\text{BIC}(\lambda) = -2\ell(\hat{\boldsymbol{\theta}}_\lambda) + k_\lambda \log(m)$$

where $m$ denotes the number of independent, declustered exceedances, and $k_\lambda$ is the number of active parameters.

---

## 2. Return Level Inference

Under non-stationarity, traditional definitions of a $T$-year return level are conceptually ill-defined. This package implements three distinct frameworks:

### 2.1 Approach A: Static Conditional Return Level
The covariates are fixed at a constant scenario $\mathbf{x}_t \equiv \mathbf{x}^{\ast}$, representing a hypothetical frozen climate state. The conditional return level $z_T(\mathbf{x}^{\ast})$ is obtained analytically:

$$z_T(\mathbf{x}^{\ast}) = \mu(\mathbf{x}^{\ast}) + \frac{\sigma(\mathbf{x}^{\ast})}{\xi(\mathbf{x}^{\ast})} \left[ \left(-\log\left(1 - \frac{1}{T}\right)\right)^{-\xi(\mathbf{x}^{\ast})} - 1 \right]$$

### 2.2 Approach B: Unconditional Parametric Stochastic Marginalization
To capture long-run risk over the natural variability of the climate system, the non-stationary intensity is integrated over the stationary joint distribution $\Pi$ of the covariate trajectories $\mathbf{v}$:

$$G_{\Pi}(z) = \mathbb{E}_{\{\mathbf{v} \sim \Pi\}} \left\lbrack G(z \mid \mathbf{v}) \right\rbrack = \mathbb{E}_{\{\mathbf{v} \sim \Pi\}} \left\lbrack \exp \left\lbrace -\frac{1}{n_{y}} \sum_{j=1}^{n_{y}} \left\lbrack 1 + \xi(t_{j}) \left\lparen \frac{z - \mu(t_{j} \mid \mathbf{v})}{\sigma(t_{j} \mid \mathbf{v})} \right\rparen \right\rbrack_{{+}}^{-1/\xi(t_{j})} \right\rbrace \right\rbrack$$

The joint distribution $\Pi$ is modeled via a stable, stationary Vector Autoregressive process, $\text{VAR}(p)$. Synthetic daily trajectories are simulated, Fourier seasonality is re-injected, and $G_{\Pi}(z)$ is estimated via Monte Carlo integration over $n_{mc}$ simulated years:

$$\hat{G}_{\Pi, B}(z) = \frac{1}{n_{mc}} \sum_{r=1}^{n_{mc}} G(z \mid \mathbf{v}^{(r)})$$

The marginalized unconditional return level $z_T^{\Pi}$ is recovered numerically as the unique root satisfying:

$$\hat{G}_{\Pi, B}(z_T^{\Pi}) - \left( 1 - \frac{1}{T} \right) = 0$$

### 2.3 Approach C: Empirical Marginalization (Non-Parametric Control)
The continuous probability space $\Pi$ is replaced by the empirical historical distribution $\hat{\Pi}$. The integrated probability is estimated as the sample mean over the fully observed historical daily multivariate trajectories of length $n_{obs}$:

$$\hat{G}_{\text{emp}, C}(z) = \frac{1}{n_{obs}} \sum_{j=1}^{n_{obs}} G(z \mid \mathbf{v}_j)$$

---

## 3. Model Diagnostics: Transformed Residuals

Goodness-of-fit is assessed using the Time-Change Theorem. The $k$-th transformed residual $Z_k$ represents the integrated intensity measure between consecutive exceedance times $t_{k-1}$ and $t_k$:

$$Z_k \approx \frac{1}{n_y} \sum_{j: t_{k-1} \le t_j < t_k} \left[ 1 + \hat{\xi}(t_j) \left( \frac{u - \hat{\mu}(t_j)}{\hat{\sigma}(t_j)} \right) \right]_{+}^{-1/\hat{\xi}(t_j)}$$

Under correct model specification, the transformed residuals are independent and identically distributed, $Z_k \sim \text{Exp}(1)$, satisfying $\mathbb{E}[Z_k] = 1$ and $\text{Var}(Z_k) = 1$.

---

## 4. Installation

You can install the stable release version of `margEVT` from CRAN:

```r
install.packages("margEVT")
```

Alternatively, you can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("rodrigo-villa/margEVT")
```

---

## 5. Quick-Start Example

This reproducible example simulates a generic non-stationary extreme value process, fits the regularized NHPP model using the LASSO penalty ($\alpha = 1$) with BIC-based $\lambda$ selection, and computes return levels under all three major marginalization frameworks (including the stochastically marginalized Approach B).

```r
library(margEVT)

# =============================================================================
# 1. Simulate Non-Stationary Extreme Value Data
# =============================================================================
set.seed(3)
n_years    <- 10L
n_per_year <- 100L
n          <- n_years * n_per_year
years      <- rep(2011:(2011 + n_years - 1L), each = n_per_year)

# Simulate 4 generic standardized covariates
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
x3 <- rnorm(n, mean = 0, sd = 1)
x4 <- rnorm(n, mean = 0, sd = 1)

# Define non-stationary GEV parameters
# Location varies with x1 and x3; Scale (log) varies with x2 and x4
mu    <- 15 + 2.0 * x1 + 1.5 * x3
sigma <- exp(1.2 + 0.3 * x2 + 0.1 * x4)
xi    <- -0.10  # Bounded upper tail (Weibull domain)

# Generate non-stationary GEV observations via the inverse transform method
u_rand <- runif(n)
y      <- mu + sigma * (((-log(u_rand))^(-xi) - 1) / xi)

sim_data <- data.frame(
  y    = y,
  x1   = x1,
  x2   = x2,
  x3   = x3,
  x4   = x4,
  year = years
)

# =============================================================================
# 2. Fit the Penalized Non-Stationary NHPP Model
# =============================================================================

fit <- fit_nhpp(
  df             = sim_data, 
  threshold      = 19,
  loc_vars       = c("x1", "x2", "x3"),
  scale_vars     = c("x1", "x2", "x4"),
  shape_vars     = NULL,       # Stationary shape parameter
  penalty        = "lasso",
  alpha          = 1.0,
  lambda         = "bic",
  obs_per_year   = 100,        # 100 observations per annual block
  calc_hessian   = TRUE,
  verbose        = TRUE
)

# Inspect estimated parameter coefficients and Hessian conditioning
summary(fit)

# Extract raw coefficients
coef(fit)

# =============================================================================
# 4. Fit the Parametric VAR(p) Stochastic Generator
# =============================================================================
# Automatically deseasonalizes active covariates, selects the optimal lag p 
# via BIC, and fits the multivariate VAR model to the anomalies.
generator <- fit_var_generator(fit, sim_data)

# Simulate 50 years of stochastic trajectories from the fitted generator,
# automatically re-injecting the seasonal harmonics.
mc_sample <- simulate_covariates(generator, n_mc = 50L, n_obs = 100L, seed = 123L)

# =============================================================================
# 5. Compute Unconditional Return Levels (Core Contribution)
# =============================================================================
# Rather than conditioning on a "frozen" covariate state (Approach A),
# we integrate out covariate uncertainty under both the parametric generator
# (Approach B) and the non-parametric empirical block bootstrap (Approach C).
rl_long <- marginalize(
  fit         = fit,
  data        = sim_data,
  TRs         = c(2, 5, 10, 50), # Return periods (years)
  n_obs       = 100L,            # 100 observations per year
  approaches  = c("A", "B", "C"),# Evaluate all three approaches
  scenarios   = list("baseline" = list(x1 = 0, x2 = 0, x3 = 0, x4 = 0)),
  mc_sample   = mc_sample,       # Supplied simulated sample from our VAR(p)
  n_boot      = 100L,            # 100 bootstrap years for Approach C
  seed        = 42L,
  year_col    = "year"
)

# Convert the long-format return levels into a wide-format risk table
rl_wide <- rl_table(rl_long)
print(rl_wide)
```

---

## 6. Citation

To cite `margEVT` in publications, please use:

> Villa, R. F. (2026). *A Novel Regularized Point Process and Stochastic Marginalization Framework for Return Level Inference under Covariate-Driven Extremes* (Master's dissertation, Instituto de Matemática e Estatística, Universidade Federal do Rio Grande do Sul, Porto Alegre, Brazil. Advisor: Flavio Ziegelmann).

Alternatively, run `citation("margEVT")` in R once the package is installed.
