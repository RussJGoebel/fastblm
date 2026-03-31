#' Apply a linear constraint to a fitted model
#'
#' Given a fastblm_fit, enforces Cx = 0 via a Schur complement correction.
#' Returns a new fastblm_fit with updated posterior mean and cached quantities
#' for constrained posterior_se.
#'
#' The constrained posterior is:
#'   mu_c     = mu - Sigma C' (C Sigma C')^{-1} C mu
#'   Sigma_c  = Sigma - Sigma C' (C Sigma C')^{-1} C Sigma
#'
#' @param fit fastblm_fit object
#' @param C q x p constraint matrix. Enforces C x = 0.
#'
#' @return fastblm_fit object with updated posterior_mean and constraint
#'   quantities cached for posterior_se.
#' @export
constrain <- function(fit, C) {
  stopifnot(inherits(fit, "fastblm_fit"))
  if (!is_matrix(C)) stop("`C` must be a matrix.")
  p <- length(fit$posterior_mean)
  if (ncol(C) != p) stop(sprintf("ncol(C) = %d but p = %d.", ncol(C), p))

  q <- nrow(C)

  # --- Step 1: compute Sigma C' = sigma2e * K^{-1} C' ------------------
  # Each column of C' is a p-vector; solve K against each one.
  # Result is p x q matrix.
  SigmaCt <- .solve_posterior(fit, t(C))   # p x q

  # --- Step 2: C Sigma C' -- q x q, small, factor it -------------------
  CSigmaCt <- C %*% SigmaCt               # q x q
  CSigmaCt <- 0.5 * (CSigmaCt + t(CSigmaCt))   # symmetrize
  CC       <- chol(as.matrix(CSigmaCt))    # upper triangular: CC' CC = CSigmaCt

  # --- Step 3: update posterior mean ------------------------------------
  # mu_c = mu - Sigma C' (C Sigma C')^{-1} C mu
  Cmu <- as.numeric(C %*% fit$posterior_mean)   # q vector
  # solve (C Sigma C') v = C mu  via CC^T CC v = Cmu
  v   <- backsolve(CC, forwardsolve(t(CC), Cmu))
  mu_c     <- fit$posterior_mean - as.numeric(SigmaCt %*% v)

  # --- Return updated fit with cached constraint quantities -------------
  fit$posterior_mean  <- mu_c
  fit$constraint      <- list(
    C        = C,
    SigmaCt  = SigmaCt,   # p x q -- needed for posterior_se correction
    CC       = CC          # Cholesky of C Sigma C' -- needed for posterior_se
  )
  fit
}

# Internal: solve K^{-1} B scaled by sigma2e, dispatching on solver type
# B is p x q matrix, returns p x q matrix Sigma B = sigma2e K^{-1} B
.solve_posterior <- function(fit, B) {
  B <- as.matrix(B)

  result <- switch(fit$solver_type,

                   cholesky = {
                     if (is.null(fit$chol_factor)) stop("No Cholesky factor in fit.")
                     as.matrix(Matrix::solve(fit$chol_factor, B)) * fit$sigma2e
                   },

                   woodbury = {
                     # Sigma = phi Q^{-1} - phi^2 Q^{-1} A' M^{-1} A Q^{-1}
                     # Apply to each column of B
                     if (is.null(fit$chol_M) || is.null(fit$QinvAt)) stop("No Woodbury quantities in fit.")
                     phi    <- fit$phi
                     QinvB  <- apply(B, 2, fit$apply_Qinv)                          # p x q
                     MinvAQinvB <- as.matrix(Matrix::solve(fit$chol_M, t(fit$QinvAt) %*% QinvB))  # n x q -- wait
                     # A Q^{-1} B is n x q
                     AQinvB     <- fit$A %*% QinvB                                   # n x q
                     MinvAQinvB <- as.matrix(Matrix::solve(fit$chol_M, AQinvB))      # n x q
                     (phi * QinvB - phi^2 * fit$QinvAt %*% MinvAQinvB) * fit$sigma2e
                   },

                   pcg = {
                     if (is.null(fit$apply_K)) stop("No apply_K in fit.")
                     apply(B, 2, function(b) pcg(fit$apply_K, b)$x) * fit$sigma2e
                   },

                   stop("Unknown solver_type.")
  )

  result
}
