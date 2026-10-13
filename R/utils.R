# =============================================================================
# utils.R
# Small shared utilities and lightweight exported helpers.
# =============================================================================

# Runs `code` with a fixed seed and restores the caller's random-number state
# afterwards (or removes the seed if none existed).
.with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

# Seeded evaluation when `seed` is supplied; otherwise uses the current stream.
.maybe_seeded <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  .with_preserved_seed(seed, code)
}

# Saves the caller's RNG kind and state; the returned function restores them.
.rng_snapshot <- function() {
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed)
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  function() {
    suppressWarnings(RNGkind(kind[1L], kind[2L], kind[3L]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

# `n` independent L'Ecuyer-CMRG streams (parallel::nextRNGStream). Replicate
# b always uses stream b, so results do not depend on how replicates are
# distributed over cores. With seed = NULL the base seed is drawn from the
# caller's stream, which therefore advances by one draw.
.rng_streams <- function(n, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  restore <- .rng_snapshot()
  on.exit(restore(), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

# Evaluates `code` on a given RNG stream and restores the caller's RNG.
.with_rng_stream <- function(stream, code) {
  restore <- .rng_snapshot()
  on.exit(restore(), add = TRUE)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  force(code)
}

.validate_seed <- function(seed, fun) {
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)))
    stop(sprintf("%s: `seed` must be NULL or a single finite number.", fun))
  invisible(seed)
}

.validate_cores <- function(cores, fun) {
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1)
    stop(sprintf("%s: `cores` must be a positive integer.", fun))
  as.integer(cores)
}

# lapply() that optionally distributes work over a PSOCK cluster. The worker
# function must be deterministic given its input, so results do not depend on
# the number of cores. Workers receive the caller's library paths, which lets
# them load margEVT from non-default libraries.
.par_lapply <- function(X, FUN, cores = 1L) {
  if (cores <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  cl <- parallel::makeCluster(min(cores, length(X)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # The setter must not be enclosed by the margEVT namespace: unserializing
  # such a closure would load margEVT on the worker before the paths are set.
  set_paths <- function(paths) .libPaths(paths)
  environment(set_paths) <- globalenv()
  parallel::clusterCall(cl, set_paths, .libPaths())
  parallel::parLapply(cl, X, FUN)
}

.penalized_parameter_indices <- function(fit) {
  dm <- fit$dm
  if (is.null(dm) || is.null(dm$X_mu) || is.null(dm$X_sigma) || is.null(dm$X_xi))
    return(NULL)

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  unique(c(
    dm$idx_pen_mu,
    p_mu + dm$idx_pen_sigma,
    if (isTRUE(fit$penalize_shape)) p_mu + p_sig + dm$idx_pen_xi else integer(0L)
  ))
}

.active_parameter_mask <- function(fit, tol) {
  pen_idx <- .penalized_parameter_indices(fit)
  if (is.null(pen_idx)) return(abs(fit$par) > tol)
  keep <- rep(TRUE, length(fit$par))
  keep[pen_idx] <- abs(fit$par[pen_idx]) > tol
  keep
}

.fit_tol <- function(fit, tol, fun) {
  tol <- if (is.null(tol)) {
    if (is.null(fit$active_tol)) 1e-2 else fit$active_tol
  } else tol
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop(sprintf("%s: `tol` must be a single positive numeric value.", fun))
  tol
}

# Names of the seasonal harmonic columns recognised by margEVT: cos<k>,
# sin<k> and the legacy spelling sen<k>, for integer k >= 1.
.harmonic_pattern <- "^(cos|sin|sen)[1-9][0-9]*$"

.is_harmonic <- function(x) grepl(.harmonic_pattern, x)


#' Compute BIC for a fitted nhpp_fit model
#'
#' Computes \eqn{-2\ell + k \log m}, where \eqn{\ell} is the unpenalized
#' log-likelihood at the operational estimator, \eqn{k} the number of active
#' parameters, and \eqn{m} the number of exceedances. Equivalent to
#' \code{BIC(fit)} for fits produced by \code{margEVT} 0.4.0 or later.
#'
#' @param fit An \code{nhpp_fit} object.
#' @param tol Numeric. Coefficients smaller than this are counted as zero
#'   for the active parameter count. If \code{NULL}, inherits
#'   \code{fit$active_tol}; legacy objects fall back to \code{1e-2}.
#'
#' @return Numeric scalar. BIC value.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)), x = rnorm(330))
#' fit <- fit_nhpp(df, threshold = 4, loc_vars = "x", penalty = "none",
#'                 verbose = FALSE)
#' bic_nhpp(fit)
#'
#' @export
bic_nhpp <- function(fit, tol = NULL) {
  if (!inherits(fit, "nhpp_fit"))
    stop("bic_nhpp: `fit` must be an nhpp_fit object.")
  nllh_bic <- if (!is.null(fit$nllh_oper)) fit$nllh_oper else fit$nllh_raw
  if (!is.finite(nllh_bic))
    return(NA_real_)
  tol <- .fit_tol(fit, tol, "bic_nhpp")
  k_active <- sum(.active_parameter_mask(fit, tol))
  2 * nllh_bic + k_active * log(fit$n_exc)
}


#' Check if an object is an nhpp_fit
#'
#' @param x Any R object.
#' @return Logical.
#' @examples
#' is_nhpp_fit(list())
#' @export
is_nhpp_fit <- function(x) inherits(x, "nhpp_fit")


#' Extract the number of exceedances from a fitted model
#'
#' @param fit An \code{nhpp_fit} object.
#' @param y Numeric vector of observations. If missing, the series stored in
#'   the fitted object is used.
#'
#' @return Integer. Number of observations exceeding the threshold.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)))
#' fit <- fit_nhpp(df, threshold = 4, penalty = "none", verbose = FALSE)
#' n_exceedances(fit)
#'
#' @export
n_exceedances <- function(fit, y) {
  if (!inherits(fit, "nhpp_fit"))
    stop("n_exceedances: `fit` must be an nhpp_fit object.")
  if (missing(y)) {
    if (is.null(fit$y))
      stop("n_exceedances: `y` must be supplied for fits without a stored series.")
    y <- fit$y
  }
  sum(y > fit$threshold, na.rm = TRUE)
}


#' Wide-format return level table from marginalize() output
#'
#' Converts the long-format output of \code{\link{marginalize}} into a
#' wide table with one row per scenario and one column per return period.
#'
#' @param marg_result Data frame returned by \code{\link{marginalize}}.
#'
#' @return A wide data frame.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = c(rexp(300), runif(30, 5, 10)),
#'                  year = rep(1991:2020, each = 11))
#' fit <- fit_nhpp(df, threshold = 4, penalty = "none", obs_per_year = 11,
#'                 verbose = FALSE)
#' rl_table(marginalize(fit, df, TRs = c(10, 50), approaches = "A"))
#'
#' @export
rl_table <- function(marg_result) {
  if (!is.data.frame(marg_result))
    stop("rl_table: `marg_result` must be a data frame from marginalize().")
  if (!all(c("approach", "scenario", "TR", "RL") %in% names(marg_result)))
    stop("rl_table: expected columns: approach, scenario, TR, RL.")

  marg_result <- as.data.frame(marg_result)
  keys     <- unique(marg_result[, c("approach", "scenario")])
  tr_vals  <- sort(unique(marg_result$TR))
  tr_cols  <- paste0("T", tr_vals)

  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    app <- keys$approach[i]
    sc  <- keys$scenario[i]
    sub <- marg_result[marg_result$approach == app &
                         marg_result$scenario == sc, ]
    rl_vals <- sub$RL[match(tr_vals, sub$TR)]
    row <- data.frame(approach = app, scenario = sc)
    for (j in seq_along(tr_cols)) row[[tr_cols[j]]] <- rl_vals[j]
    row
  }))
}
