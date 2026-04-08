#' Project out fixed effects for REML
#'
#' Projects \eqn{y} and \eqn{A} onto the orthogonal complement of
#' \eqn{X_{\mathrm{fixed}}} with respect to the \eqn{R^{-1}}-weighted inner
#' product. When \code{R_inv} is NULL (identity), this reduces to ordinary QR
#' residualisation.
#'
#' For general \eqn{R}, the projection whitens by \eqn{R^{-1/2}},
#' residualises, then back-transforms, so the returned \eqn{y^*} and
#' \eqn{A^*} satisfy \eqn{X^\top R^{-1} y^* = 0} and
#' \eqn{X^\top R^{-1} A^* = 0}.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param X_fixed n x q matrix of fixed covariates
#' @param R_inv n x n inverse noise covariance matrix. NULL = identity.
#'
#' @return list with projected y, A, and n_eff
#' @export
reml_project <- function(y, A, X_fixed, R_inv = NULL) {
  if (!is_matrix(X_fixed)) stop("`X_fixed` must be a matrix.")
  n <- length(y)
  if (nrow(X_fixed) != n) stop("nrow(X_fixed) must equal length(y).")

  if (is.null(R_inv)) {
    # R = I: plain QR residualisation
    qr_fixed <- qr(as.matrix(X_fixed))
    q_rank   <- qr_fixed$rank
    if (q_rank < ncol(X_fixed))
      warning(sprintf("X_fixed has rank %d < ncol = %d.", q_rank, ncol(X_fixed)))
    return(list(
      y     = qr.resid(qr_fixed, y),
      A     = qr.resid(qr_fixed, as.matrix(A)),
      n_eff = n - q_rank
    ))
  }

  # General R: whiten, project, back-transform.
  # R^{-1} = L L' (Cholesky), so R^{-1/2} = L.
  # Whiten: y_w = L y,  A_w = L A,  X_w = L X.
  # Residualise in whitened space, then return un-whitened residuals.
  Rinv_mat <- as.matrix(R_inv)
  L <- t(chol(Rinv_mat))           # lower-triangular, L L' = R^{-1}
  y_w <- as.numeric(L %*% y)
  A_w <- L %*% as.matrix(A)
  X_w <- L %*% as.matrix(X_fixed)

  qr_fixed <- qr(X_w)
  q_rank   <- qr_fixed$rank
  if (q_rank < ncol(X_fixed))
    warning(sprintf("X_fixed has rank %d < ncol = %d.", q_rank, ncol(X_fixed)))

  # Residuals in whitened space, then back-transform by L^{-1}
  y_w_res <- qr.resid(qr_fixed, y_w)
  A_w_res <- qr.resid(qr_fixed, A_w)
  Linv    <- solve(L)
  list(
    y     = as.numeric(Linv %*% y_w_res),
    A     = Linv %*% A_w_res,
    n_eff = n - q_rank
  )
}
