# =============================================================================
# decluster.R
# Extremal index and declustering of threshold exceedances.
#
#   - extremal_index() : intervals (Ferro and Segers, 2003) or runs estimator
#   - decluster()      : keeps the largest value of each cluster
#
# Clusters are separated by more than `run_length` consecutive
# non-exceedances. For the intervals method the run length is the one
# implied by the intervals estimate (Ferro and Segers, 2003, Sec. 5).
# =============================================================================

# Exceedance times and inter-exceedance times of a series.
.exceedance_times <- function(y, threshold, fun) {
  if (!is.numeric(y) || length(y) < 2L || any(!is.finite(y)))
    stop(sprintf("%s: `y` must be a numeric vector of finite values.", fun))
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold))
    stop(sprintf("%s: `threshold` must be a single finite value.", fun))
  s <- which(y > threshold)
  if (length(s) < 2L)
    stop(sprintf("%s: fewer than two exceedances of `threshold`.", fun))
  list(s = s, gaps = diff(s))
}

# Intervals estimator of the extremal index (Ferro and Segers, 2003, Eq. 4):
# theta_hat when every gap is at most 2, the bias-corrected theta_tilde
# otherwise; truncated at 1.
.theta_intervals <- function(gaps) {
  N <- length(gaps) + 1L
  theta <- if (max(gaps) <= 2)
    2 * sum(gaps)^2 / ((N - 1) * sum(gaps^2)) else
      2 * sum(gaps - 1)^2 / ((N - 1) * sum((gaps - 1) * (gaps - 2)))
  min(1, theta)
}

# Run length implied by an extremal index estimate: the C - 1 largest gaps
# separate clusters, C = ceiling(theta * N), so the run length is the Cth
# largest gap and clusters are split by the gaps above it. Gaps tied with
# the Cth largest do not split, so ties can give fewer than C clusters.
.run_length_intervals <- function(gaps, theta) {
  N <- length(gaps) + 1L
  C <- max(1L, min(N, as.integer(ceiling(theta * N))))
  if (C >= N) return(0L)
  as.integer(sort(gaps, decreasing = TRUE)[C])
}

# Cluster labels of the exceedances for a given run length.
.cluster_labels <- function(gaps, run_length)
  c(1L, 1L + cumsum(gaps > run_length))

#' Extremal index of a series of threshold exceedances
#'
#' Estimates the extremal index \eqn{\theta}, the reciprocal of the mean
#' cluster size of threshold exceedances (\eqn{\theta = 1} for no
#' clustering).
#'
#' \describe{
#'   \item{\code{"intervals"}}{The intervals estimator of Ferro and Segers
#'     (2003), based on the first two moments of the inter-exceedance times
#'     \eqn{T_i}: \eqn{2(\sum T_i)^2 / \{(N-1)\sum T_i^2\}} when
#'     \eqn{\max T_i \le 2}, and
#'     \eqn{2\{\sum (T_i - 1)\}^2 / \{(N-1)\sum (T_i-1)(T_i-2)\}} otherwise,
#'     truncated at 1. It needs no run length.}
#'   \item{\code{"runs"}}{The number of clusters over the number of
#'     exceedances, clusters being separated by more than \code{run_length}
#'     consecutive non-exceedances.}
#' }
#'
#' @param y Numeric vector, the series in time order (regularly spaced, with
#'   no missing values).
#' @param threshold Single numeric threshold.
#' @param method \code{"intervals"} (default) or \code{"runs"}.
#' @param run_length Positive integer, the run length for
#'   \code{method = "runs"}. Default \code{1L}.
#'
#' @return A named numeric vector with \code{theta}, \code{n_clusters} and
#'   \code{run_length}. For the intervals method these are the clusters and
#'   run length used by \code{\link{decluster}}.
#'
#' @references
#' Ferro, C. A. T. and Segers, J. (2003). Inference for clusters of extreme
#' values. \emph{Journal of the Royal Statistical Society: Series B}, 65(2),
#' 545--556.
#'
#' @seealso \code{\link{decluster}}
#'
#' @examples
#' set.seed(1)
#' y <- as.numeric(stats::arima.sim(list(ar = 0.7), 5000))
#' u <- stats::quantile(y, 0.95)
#' extremal_index(y, u)
#' extremal_index(y, u, method = "runs", run_length = 3)
#'
#' @export
extremal_index <- function(y, threshold, method = c("intervals", "runs"),
                           run_length = 1L) {
  method <- match.arg(method)
  et <- .exceedance_times(y, threshold, "extremal_index")
  N  <- length(et$s)
  if (method == "intervals") {
    theta <- .theta_intervals(et$gaps)
    r     <- .run_length_intervals(et$gaps, theta)
  } else {
    if (!is.numeric(run_length) || length(run_length) != 1L ||
        !is.finite(run_length) || run_length < 0)
      stop("extremal_index: `run_length` must be a non-negative integer.")
    r <- as.integer(run_length)
  }
  n_cl <- max(.cluster_labels(et$gaps, r))
  if (method == "runs") theta <- n_cl / N
  c(theta = theta, n_clusters = n_cl, run_length = r)
}

#' Decluster a series of threshold exceedances
#'
#' Groups the exceedances of \code{threshold} into clusters separated by more
#' than \code{run_length} consecutive non-exceedances, keeps the largest
#' value of each cluster (the first, in case of ties) and replaces the other
#' exceedances by \code{replace_with}, by default the threshold itself, so
#' that they no longer exceed it. Values below the threshold are unchanged,
#' and the result can be passed directly to \code{\link{fit_nhpp}} or used as
#' the \code{decluster} function of \code{\link{backtest}}.
#'
#' With \code{method = "intervals"} (default), the run length is estimated
#' from the data by the intervals method of Ferro and Segers (2003): the
#' \eqn{C - 1} largest inter-exceedance times separate clusters, with
#' \eqn{C = \lceil \hat\theta N \rceil}. With \code{"runs"}, the run length
#' is \code{run_length}.
#'
#' @inheritParams extremal_index
#' @param run_length Positive integer, the run length for
#'   \code{method = "runs"}. Ignored for \code{"intervals"}.
#' @param replace_with Value given to the exceedances that are not cluster
#'   maxima. Default: \code{threshold}.
#'
#' @return The declustered series, with attributes \code{theta},
#'   \code{run_length}, \code{n_clusters} and \code{clusters} (the cluster
#'   label of each exceedance of the original series).
#'
#' @note \code{extRemes::decluster(method = "intervals")} (version 2.2-1)
#'   estimates the run length but then declusters with a run length of 1;
#'   \code{decluster(method = "runs", run_length = 1)} reproduces that output.
#'
#' @references
#' Ferro, C. A. T. and Segers, J. (2003). Inference for clusters of extreme
#' values. \emph{Journal of the Royal Statistical Society: Series B}, 65(2),
#' 545--556.
#'
#' @seealso \code{\link{extremal_index}}
#'
#' @examples
#' set.seed(1)
#' y <- as.numeric(stats::arima.sim(list(ar = 0.7), 5000))
#' u <- stats::quantile(y, 0.95)
#' yd <- decluster(y, u)
#' c(exceedances = sum(y > u), clusters = sum(yd > u))
#' attr(yd, "run_length")
#'
#' @export
decluster <- function(y, threshold, method = c("intervals", "runs"),
                      run_length = 1L, replace_with = threshold) {
  method <- match.arg(method)
  ei <- extremal_index(y, threshold, method = method, run_length = run_length)
  if (!is.numeric(replace_with) || length(replace_with) != 1L ||
      is.na(replace_with) || replace_with > threshold)
    stop("decluster: `replace_with` must be a single value not above `threshold`.")
  et <- .exceedance_times(y, threshold, "decluster")
  cl <- .cluster_labels(et$gaps, ei[["run_length"]])
  # first position of the maximum within each cluster
  vals <- y[et$s]
  is_max <- stats::ave(vals, cl, FUN = function(v) seq_along(v) == which.max(v))
  out <- as.numeric(y)
  out[et$s[is_max == 0]] <- replace_with
  attr(out, "theta")      <- ei[["theta"]]
  attr(out, "run_length") <- ei[["run_length"]]
  attr(out, "n_clusters") <- ei[["n_clusters"]]
  attr(out, "clusters")   <- cl
  out
}
