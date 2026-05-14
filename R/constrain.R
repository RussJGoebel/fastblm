#' Apply a linear constraint to a fitted model
#'
#' Given a fastblm_fit, enforces \eqn{Cx = 0} via a Schur complement
#' correction. Returns a new fastblm_fit with updated posterior mean and
#' cached quantities for constrained \code{posterior_se}.
#'
#' The constrained posterior is:
#' \deqn{\mu_c = \mu - \Sigma C^\top (C \Sigma C^\top)^{-1} C \mu}
#' \deqn{\Sigma_c = \Sigma - \Sigma C^\top (C \Sigma C^\top)^{-1} C \Sigma}
#'
#' Numerically, the mean correction is computed via \eqn{K^{-1} C^\top}
#' (without \eqn{\sigma^{2}_{e}}) since \eqn{\sigma^{2}_{e}} cancels exactly in
#' \eqn{\Sigma C^\top (C \Sigma C^\top)^{-1}}. This avoids catastrophic
#' cancellation when \eqn{\sigma^{2}_{e}} is small and \eqn{K^{-1}} is large.
#' \eqn{\sigma^{2}_{e}} is reintroduced only for the cached \code{SigmaCt} used
#' by \code{posterior_se}.
#'
#' @param fit fastblm_fit object
#' @param C q x p constraint matrix. Enforces \eqn{Cx = 0}.
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

  # --- Step 1: compute K^{-1} C' (WITHOUT sigma2e) -------------------------
  # sigma2e cancels in the mean correction:
  #   Sigma C' (C Sigma C')^{-1} = sigma2e K^{-1} C' (C sigma2e K^{-1} C')^{-1}
  #                               = K^{-1} C' (C K^{-1} C')^{-1}
  # Working without sigma2e avoids catastrophic cancellation when sigma2e is
  # tiny (~1e-6) and K^{-1} is large (intrinsic prior, large phi).
  KinvCt <- .solve_posterior(fit, t(C)) / fit$sigma2e   # p x q, pure K^{-1} C'

  # --- Step 2: C K^{-1} C' -- q x q, factor it -----------------------------
  CKinvCt <- C %*% KinvCt
  CKinvCt <- 0.5 * (CKinvCt + t(CKinvCt))              # symmetrize
  CC_Kinv <- chol(as.matrix(CKinvCt))                   # upper triangular

  # --- Step 3: update posterior mean (sigma2e-free, exact cancellation) -----
  # mu_c = mu - K^{-1} C' (C K^{-1} C')^{-1} C mu
  Cmu  <- as.numeric(C %*% fit$posterior_mean)           # q vector
  v    <- backsolve(CC_Kinv, forwardsolve(t(CC_Kinv), Cmu))
  mu_c <- fit$posterior_mean - as.numeric(KinvCt %*% v)

  # --- Step 4: cache Sigma C' and Chol(C Sigma C') for posterior_se ---------
  # posterior_se needs sigma2e in SigmaCt to compute diag(Sigma_c) correctly.
  SigmaCt  <- KinvCt * fit$sigma2e                      # p x q = sigma2e K^{-1} C'
  CSigmaCt <- C %*% SigmaCt
  CSigmaCt <- 0.5 * (CSigmaCt + t(CSigmaCt))
  CC       <- chol(as.matrix(CSigmaCt))                 # Cholesky of C Sigma C'

  # --- Return updated fit ---------------------------------------------------
  fit$posterior_mean <- mu_c
  fit$constraint     <- list(
    C       = C,
    SigmaCt = SigmaCt,  # p x q, sigma2e * K^{-1} C' -- for posterior_se
    CC      = CC         # Cholesky of C Sigma C'      -- for posterior_se
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
