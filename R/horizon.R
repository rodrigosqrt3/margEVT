# =============================================================================
# horizon.R
# Risk over a planning horizon of several years.
#
# A trajectory is a sequence of L annual covariate blocks (years 1, ..., L of
# the horizon); an ensemble is a list of S such trajectories. For year l of
# trajectory s, G_ls(z) is the annual non-exceedance probability, and the
# probability of no exceedance of z over the horizon is prod_l G_ls(z),
# averaged over the ensemble. Non-stationary scenarios (e.g. covariates from
# climate projections, or a trend running past the record) enter as
# trajectories whose years differ.
#
#   - horizon_risk()      : exceedance probability, expected number of years
#                           with an exceedance, and waiting time over L years
#   - design_life_level() : level exceeded with probability p over L years
# =============================================================================

# Validates `trajectories` (one trajectory or an ensemble of equal length)
# and returns the stacked model parameters, year-major within trajectories.
.horizon_params <- function(fit, trajectories, n_obs, period, interactions,
                            fixed, fun) {
  if (!inherits(fit, "nhpp_fit"))
    stop(sprintf("%s: `fit` must be an nhpp_fit object.", fun))
  if (.is_block_list(trajectories)) trajectories <- list(trajectories)
  if (!is.list(trajectories) || length(trajectories) < 1L ||
      !all(vapply(trajectories, .is_block_list, logical(1L))))
    stop(sprintf(paste0("%s: `trajectories` must be a list of annual data ",
                        "frames or a list of such lists."), fun))
  L <- unique(vapply(trajectories, length, integer(1L)))
  if (length(L) != 1L)
    stop(sprintf("%s: all trajectories must have the same number of years.",
                 fun))
  fixed <- .validate_fixed(fixed, fun)
  n_obs <- if (is.null(n_obs)) as.integer(round(fit$obs_per_year)) else
    as.integer(n_obs)
  if (!is.finite(n_obs) || n_obs < 1L)
    stop(sprintf("%s: `n_obs` must be a positive integer.", fun))

  cb <- .covariate_builder(fit, n_obs, period, interactions, fixed)
  blocks <- unlist(trajectories, recursive = FALSE)
  .check_blocks(blocks, cb$required, "trajectories (year)", fun)
  big <- cb$build_stack(blocks)
  list(params = predict_params(fit, big), n_obs = n_obs, L = L,
       S = length(trajectories))
}

#' Exceedance risk over a planning horizon
#'
#' For each level \eqn{z}, computes the probability that \eqn{z} is exceeded
#' at least once in a horizon of \eqn{L} years, the expected number of years
#' with an exceedance, and the expected waiting time to the first
#' exceedance, under covariate trajectories supplied by the user.
#'
#' Year \eqn{l} of trajectory \eqn{s} has annual non-exceedance probability
#' \eqn{G_{ls}(z) = \exp\{-n_y^{-1}\sum_t \tau_t(z)\}}, as in
#' \code{\link{annual_prob}}. With \eqn{S} trajectories,
#' \deqn{P(\textrm{no exceedance in } L \textrm{ years}) =
#'   S^{-1}\sum_s \prod_{l=1}^{L} G_{ls}(z).}
#' The years of a trajectory are conditionally independent given the
#' covariates, so the dependence between years enters only through the
#' covariates. For a stationary climate with independent years, the
#' probability of at least one exceedance of the \eqn{T}-year marginal level
#' in \eqn{L} years is \eqn{1 - (1 - 1/T)^L}; non-stationary scenarios (for
#' example, a trend that keeps growing or covariates from climate projections)
#' enter as trajectories whose years differ.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param z Numeric vector of levels.
#' @param trajectories Either one trajectory, a list of \eqn{L} data frames
#'   (the covariates of years \eqn{1, \ldots, L}, \code{n_obs} rows each), or
#'   an ensemble, a list of such lists of the same length (e.g. simulated
#'   scenarios). Each data frame must supply the active covariates of
#'   \code{fit}, as \code{mc_sample} in \code{\link{marginalize}}.
#' @param n_obs Integer observations per year. If \code{NULL}, the rounded
#'   \code{fit$obs_per_year}.
#' @param period Seasonal period; see \code{\link{build_cov_annual}}.
#' @param interactions,fixed As in \code{\link{marginalize}}.
#'
#' @return An object of class \code{nhpp_horizon} (a data frame, with
#'   \code{print} and \code{plot} methods; see
#'   \code{\link{nhpp_horizon-methods}}) with one row per level and columns \code{z},
#'   \code{prob_exceed} (probability of at least one exceedance in the
#'   horizon), \code{se_mc} (its Monte Carlo standard error over the
#'   ensemble; \code{NA} for a single trajectory), \code{expected_years}
#'   (expected number of years with an exceedance) and \code{wait}
#'   (\eqn{E[\min(W, L)]}, the expected number of years until the first
#'   exceedance \eqn{W}, counting at most \eqn{L}). The \eqn{L \times}
#'   \code{length(z)} matrix of probabilities that the first exceedance occurs
#'   by year \eqn{l} is attached as attribute \code{"cumulative"}.
#'
#' @seealso \code{\link{design_life_level}}, \code{\link{marginalize}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500), year = rep(1:30, each = 50))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' # 20 future years in which x rises linearly from 0 to 2
#' future <- lapply(seq(0, 2, length.out = 20),
#'                  function(v) data.frame(x = rep(v, 50)))
#' horizon_risk(fit, z = c(25, 30), trajectories = future)
#'
#' @export
horizon_risk <- function(fit, z, trajectories, n_obs = NULL, period = NULL,
                         interactions = list(), fixed = list()) {
  if (!is.numeric(z) || length(z) < 1L || any(!is.finite(z)))
    stop("horizon_risk: `z` must contain finite numeric values.")
  hp <- .horizon_params(fit, trajectories, n_obs, period, interactions,
                        fixed, "horizon_risk")
  prm <- hp$params
  cumulative <- matrix(NA_real_, hp$L, length(z),
                       dimnames = list(paste0("year", seq_len(hp$L)),
                                       format(z)))
  rows <- lapply(seq_along(z), function(i) {
    G <- matrix(.pathwise_prob(z[i], prm$mu, prm$sigma, prm$xi, hp$n_obs),
                nrow = hp$L)
    surv <- apply(G, 2L, cumprod)
    surv <- matrix(surv, nrow = hp$L)
    none <- surv[hp$L, ]
    cumulative[, i] <<- 1 - rowMeans(surv)
    data.frame(
      z              = z[i],
      prob_exceed    = 1 - mean(none),
      se_mc          = if (hp$S > 1L) stats::sd(none) / sqrt(hp$S) else
        NA_real_,
      expected_years = mean(colSums(1 - G)),
      wait           = 1 + sum(rowMeans(surv)[-hp$L])
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "cumulative") <- cumulative
  attr(out, "horizon") <- hp$L
  attr(out, "n_trajectories") <- hp$S
  class(out) <- c("nhpp_horizon", class(out))
  out
}

#' Design life level over a planning horizon
#'
#' The design life level (Rootzen and Katz, 2013) is the level whose
#' probability of being exceeded at least once during a horizon of \eqn{L}
#' years is \code{prob}: it solves
#' \deqn{S^{-1}\sum_s \prod_{l=1}^{L} G_{ls}(z) = 1 - p}
#' over covariate trajectories (see \code{\link{horizon_risk}}). For
#' \eqn{L = 1} and a stationary covariate law it is the marginal return level
#' of period \eqn{1/p} computed by \code{\link{marginalize}}.
#'
#' @inheritParams horizon_risk
#' @param prob Numeric vector of exceedance probabilities over the horizon,
#'   in (0, 1). Default \code{c(0.5, 0.1, 0.05)}.
#' @param tol_root Root-finding tolerance. Default \code{1e-4}.
#' @param z_hi Upper search bound. If \code{NULL},
#'   \code{3 * max(fit$y)}.
#'
#' @return A data frame with columns \code{prob}, \code{horizon} (\eqn{L}),
#'   \code{DLL} and \code{se_mc}, the Monte Carlo standard error of the root
#'   over the ensemble (\code{NA} for a single trajectory).
#'
#' @references
#' Rootzen, H. and Katz, R. W. (2013). Design life level: quantifying risk in
#' a changing climate. \emph{Water Resources Research}, 49(9), 5964--5972.
#'
#' @seealso \code{\link{horizon_risk}}, \code{\link{marginalize}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(1500), year = rep(1:30, each = 50))
#' df$y <- 10 + df$x + rexp(1500, 1 / 3)
#' fit <- fit_nhpp(df, threshold = 18, loc_vars = "x", penalty = "none",
#'                 obs_per_year = 50, verbose = FALSE)
#' # ensemble of 200 horizons of 10 years, resampling the observed years
#' traj <- replicate(200, lapply(sample(1:30, 10, replace = TRUE),
#'                               function(yr) df[df$year == yr, "x", drop = FALSE]),
#'                   simplify = FALSE)
#' design_life_level(fit, traj, prob = c(0.5, 0.1))
#'
#' @export
design_life_level <- function(fit, trajectories, prob = c(0.5, 0.1, 0.05),
                              n_obs = NULL, period = NULL,
                              interactions = list(), fixed = list(),
                              tol_root = 1e-4, z_hi = NULL) {
  if (!is.numeric(prob) || length(prob) < 1L || any(!is.finite(prob)) ||
      any(prob <= 0) || any(prob >= 1))
    stop("design_life_level: `prob` must lie strictly between 0 and 1.")
  hp <- .horizon_params(fit, trajectories, n_obs, period, interactions,
                        fixed, "design_life_level")
  z_lo <- fit$threshold
  if (is.null(z_hi)) {
    if (is.null(fit$y))
      stop("design_life_level: supply `z_hi` for fits without a stored series.")
    z_hi <- 3 * max(fit$y)
  }
  if (!is.numeric(z_hi) || length(z_hi) != 1L || !is.finite(z_hi) ||
      z_hi <= z_lo)
    stop("design_life_level: `z_hi` must be a finite value above the threshold.")
  curve <- .level_curve(hp$params, hp$n_obs, n_rows = hp$n_obs * hp$L)
  dll <- .solve_return_levels(1 / prob, curve, z_lo, z_hi, tol_root)
  se  <- vapply(dll, curve$se, numeric(1L))
  data.frame(prob = prob, horizon = hp$L, DLL = round(dll, 4L), se_mc = se)
}
