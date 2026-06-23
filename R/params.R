# =============================================================================
# params.R
# Extract time-varying parameters from a fitted nhpp_fit object.
#
# predict_params() is the single bridge between a fitted model and all
# downstream computations (marginalization, return levels, bootstrap).
# It re-evaluates the design matrices on arbitrary new data, so it works
# for in-sample fitted values, out-of-sample prediction, and simulation.
# =============================================================================

#' Compute time-varying GEV parameters from a fitted model
#'
#' Given an \code{nhpp_fit} object and a data frame of covariate values,
#' returns the implied location, scale, and shape parameters at each row.
#'
#' @param fit An object of class \code{nhpp_fit}.
#' @param newdata A data frame with one row per time point. Must contain all
#'   columns that appear in the model's design matrices (except the intercept).
#'   If \code{NULL}, returns the fitted values stored in \code{fit$fitted}.
#'
#' @return A list with three numeric vectors of length \code{nrow(newdata)}:
#'   \item{mu}{Location parameter.}
#'   \item{sigma}{Scale parameter (always positive).}
#'   \item{xi}{Shape parameter.}
#'
#' @export
predict_params <- function(fit, newdata = NULL) {

  if (!inherits(fit, "nhpp_fit"))
    stop("predict_params: `fit` must be an nhpp_fit object.")

  # Return stored fitted values for in-sample case
  if (is.null(newdata)) {
    return(fit$fitted)
  }

  if (!is.data.frame(newdata))
    stop("predict_params: `newdata` must be a data frame.")

  dm  <- fit$dm
  par <- fit$par

  p_mu  <- ncol(dm$X_mu)
  p_sig <- ncol(dm$X_sigma)
  p_xi  <- ncol(dm$X_xi)

  beta_mu    <- par[seq_len(p_mu)]
  beta_sigma <- par[p_mu + seq_len(p_sig)]
  beta_xi    <- par[p_mu + p_sig + seq_len(p_xi)]

  # Re-evaluate design matrices on newdata using the same column structure
  # as the original fit — this is what makes out-of-sample prediction work.
  X_mu_new  <- .make_X_pred(dm$X_mu,    newdata)
  X_sig_new <- .make_X_pred(dm$X_sigma, newdata)
  X_xi_new  <- .make_X_pred(dm$X_xi,   newdata)

  list(
    mu    = as.numeric(X_mu_new  %*% beta_mu),
    sigma = exp(as.numeric(X_sig_new %*% beta_sigma)),
    xi    = as.numeric(X_xi_new  %*% beta_xi)
  )
}

# -----------------------------------------------------------------------------
# .make_X_pred() — internal
# Rebuilds a design matrix on new data using the column structure of the
# original design matrix. Handles the intercept automatically.
# -----------------------------------------------------------------------------
.make_X_pred <- function(X_template, newdata) {

  cn      <- colnames(X_template)
  cn_data <- cn[cn != "(Intercept)"]
  n       <- nrow(newdata)

  # Intercept-only model
  if (length(cn_data) == 0L) {
    return(matrix(1, nrow = n, ncol = 1L,
                  dimnames = list(NULL, "(Intercept)")))
  }

  missing_c <- setdiff(cn_data, names(newdata))
  if (length(missing_c) > 0L)
    stop(sprintf(
      "predict_params: columns required by model not found in newdata: %s",
      paste(missing_c, collapse = ", ")
    ))

  rhs <- paste(cn_data, collapse = " + ")
  stats::model.matrix(stats::as.formula(paste("~", rhs)), data = newdata)
}
