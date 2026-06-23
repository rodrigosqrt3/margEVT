# =============================================================================
# design.R
# Build design matrices for the NHPP model.
# =============================================================================

#' Build design matrices for the NHPP point process model
#'
#' @param df A data frame containing all required columns.
#' @param loc_vars Character vector of covariate names for the location (mu) parameter.
#' @param scale_vars Character vector of covariate names for the scale (sigma) parameter.
#' @param shape_vars Character vector of covariate names for the shape (xi) parameter.
#' @param free_vars Character vector of column names that should never be penalized
#'   (e.g. seasonality terms like \code{c("cos1", "sen1")}). Defaults to \code{NULL}.
#'
#' @return A list with elements:
#'   \item{X_mu}{Design matrix for location.}
#'   \item{X_sigma}{Design matrix for scale.}
#'   \item{X_xi}{Design matrix for shape.}
#'   \item{idx_pen_mu}{Integer indices of penalized columns in X_mu.}
#'   \item{idx_pen_sigma}{Integer indices of penalized columns in X_sigma.}
#'   \item{idx_pen_xi}{Integer indices of penalized columns in X_xi.}
#'
#' @export
build_design_matrices <- function(df, loc_vars = NULL, scale_vars = NULL,
                                  shape_vars = NULL, free_vars = NULL) {

  n <- nrow(df)
  never_pen <- c("(Intercept)", free_vars)

  make_X <- function(vars, label) {
    if (is.null(vars) || length(vars) == 0L) {
      return(matrix(1, nrow = n, ncol = 1L,
                    dimnames = list(NULL, "(Intercept)")))
    }
    missing_v <- setdiff(vars, names(df))
    if (length(missing_v) > 0L)
      stop(sprintf(
        "build_design_matrices [%s]: columns not found in df: %s",
        label, paste(missing_v, collapse = ", ")
      ))
    X <- as.matrix(cbind(1, df[, vars, drop = FALSE]))
    colnames(X)[1L] <- "(Intercept)"
    X
  }

  pen_idx <- function(X) {
    which(!colnames(X) %in% never_pen)
  }

  X_mu    <- make_X(loc_vars,   "mu")
  X_sigma <- make_X(scale_vars, "sigma")
  X_xi    <- make_X(shape_vars, "xi")

  list(
    X_mu          = X_mu,
    X_sigma       = X_sigma,
    X_xi          = X_xi,
    idx_pen_mu    = pen_idx(X_mu),
    idx_pen_sigma = pen_idx(X_sigma),
    idx_pen_xi    = pen_idx(X_xi)
  )
}
