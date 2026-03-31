#' REML projection
#'
#' Projects y and A onto the orthogonal complement of X_fixed,
#' removing the influence of fixed effects on hyperparameter estimation.
#'
#' @param y numeric response vector
#' @param A n x p design matrix
#' @param X_fixed n x q matrix of fixed covariates
#'
#' @return list with projected y, A, and n_eff
#' @export
reml_project <- function(y, A, X_fixed) {
  if (!is_matrix(X_fixed)) stop("X_fixed must be a matrix.")
  if (nrow(X_fixed) != length(y)) stop("nrow(X_fixed) must equal length(y).")

  qr_fixed <- qr(as.matrix(X_fixed))
  q_rank   <- qr_fixed$rank

  if (q_rank < ncol(X_fixed))
    warning(sprintf("X_fixed has rank %d < ncol = %d.", q_rank, ncol(X_fixed)))

  list(
    y     = qr.resid(qr_fixed, y),
    A     = qr.resid(qr_fixed, as.matrix(A)),
    n_eff = length(y) - q_rank
  )
}
