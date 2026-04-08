#' Posterior standard errors for coefficients or linear combinations
#'
#' Computes \eqn{\sqrt{\mathrm{diag}(A_{\mathrm{new}} \Sigma_{\mathrm{post}} A_{\mathrm{new}}^\top)}}
#' where \eqn{\Sigma_{\mathrm{post}} = \sigma^2_e K^{-1}}.
#' If \code{A_new} is NULL, returns marginal SEs for the coefficients themselves.
#' If the fit has an active constraint (from \code{constrain()}), the Schur
#' correction is applied automatically.
#'
#' Three paths depending on fit solver:
#' \describe{
#'   \item{\code{"cholesky"}}{Exact, via p x p triangular solves.}
#'   \item{\code{"woodbury"}}{Exact, via cached n x n factor. Cheap when p >> n.}
#'   \item{\code{"pcg"}}{Stochastic Hutchinson estimator, \code{n_probes} PCG solves regardless of p.}
#' }
#'
#' @param fit fastblm_fit object
#' @param A_new optional matrix of linear combinations (n_new x p). NULL = identity.
#'   May be a sparse Matrix (e.g. \code{[I_p | X_grid]} for augmented fits).
#' @param n_probes number of Hutchinson probes (PCG path only)
#'
#' @return numeric vector of posterior SEs
#' @export
posterior_se <- function(fit, A_new = NULL, n_probes = 50L) {
  stopifnot(inherits(fit, "fastblm_fit"))
  p <- length(fit$posterior_mean)

  # unconstrained diagonal variance
  diag_var <- switch(fit$solver_type,
                     cholesky = .diag_var_cholesky(fit, A_new, p),
                     woodbury = .diag_var_woodbury(fit, A_new, p),
                     pcg      = .diag_var_hutchinson(fit, A_new, p, n_probes),
                     stop("Unknown solver_type in fit object.")
  )

  # apply Schur correction if constrained
  if (!is.null(fit$constraint)) {
    diag_var <- diag_var - .constraint_correction(fit, A_new, p)
  }

  sqrt(as.numeric(diag_var))
}

# -----------------------------------------------------------------------------
# Diagonal variance paths (return diag_var, not SE)
# -----------------------------------------------------------------------------

# Cholesky path: exact
# diag(A K^{-1} A') = colSums(A' * (K^{-1} A'))
.diag_var_cholesky <- function(fit, A_new, p) {
  C <- fit$chol_factor
  if (is.null(C)) stop("No Cholesky factor found in fit object.")
  A <- if (is.null(A_new)) Matrix::Diagonal(p) else A_new
  Z <- Matrix::solve(C, Matrix::t(A))
  as.numeric(Matrix::colSums(Matrix::t(A) * Z)) * fit$sigma2e
}

# Woodbury path: exact
# Sigma = phi Q^{-1} - phi^2 Q^{-1} A' M^{-1} A Q^{-1}
.diag_var_woodbury <- function(fit, A_new, p) {
  CM     <- fit$chol_M
  QinvAt <- fit$QinvAt
  phi    <- fit$phi
  if (is.null(CM) || is.null(QinvAt))
    stop("Woodbury cached quantities not found in fit object.")

  MinvAQinv <- as.matrix(Matrix::solve(CM, t(QinvAt)))   # n x p

  if (is.null(A_new)) {
    diag_Qinv <- vapply(seq_len(p), function(i) {
      ei <- rep(0, p); ei[i] <- 1
      fit$apply_Qinv(ei)[i]
    }, numeric(1L))
    diag_corr <- Matrix::rowSums(QinvAt * t(MinvAQinv))
    (phi * diag_Qinv - phi^2 * diag_corr) * fit$sigma2e

  } else {
    QinvAtnew    <- apply(t(A_new), 2, fit$apply_Qinv)
    diag_AQinvAt <- Matrix::colSums(t(A_new) * QinvAtnew)
    AnewQinvAt   <- A_new %*% QinvAt
    MinvAnew     <- as.matrix(Matrix::solve(CM, t(AnewQinvAt)))
    diag_corr    <- Matrix::rowSums(AnewQinvAt * t(MinvAnew))
    (phi * diag_AQinvAt - phi^2 * diag_corr) * fit$sigma2e
  }
}

# PCG/Hutchinson path: stochastic
# n_probes solves regardless of p or n_new
.diag_var_hutchinson <- function(fit, A_new, p, n_probes) {
  if (is.null(fit$apply_K)) stop("No apply_K found in fit object.")
  apply_Kinv <- function(v) pcg(fit$apply_K, v)$x

  if (is.null(A_new)) {
    probes <- .rademacher(p, n_probes)
    .hutchinson_diag(apply_Kinv, probes) * fit$sigma2e
  } else {
    A     <- A_new
    n_new <- nrow(A)
    probes <- .rademacher(n_new, n_probes)
    apply_AKinvAt <- function(z) {
      as.numeric(A %*% apply_Kinv(as.numeric(Matrix::crossprod(A, z))))
    }
    .hutchinson_diag(apply_AKinvAt, probes) * fit$sigma2e
  }
}

# -----------------------------------------------------------------------------
# Constraint Schur correction
# diag(Sigma C' (C Sigma C')^{-1} C Sigma)
# = rowSums((A_new SigmaCt L^{-T})^2)  where L L' = C Sigma C'
# -----------------------------------------------------------------------------
.constraint_correction <- function(fit, A_new, p) {
  SigmaCt <- fit$constraint$SigmaCt   # p x q
  CC      <- fit$constraint$CC        # upper Cholesky of C Sigma C'

  # W = SigmaCt (CC^T)^{-1}  -- p x q
  W <- t(forwardsolve(t(CC), t(SigmaCt)))   # p x q  (dense, q is small)

  if (is.null(A_new)) {
    # W is dense p x q -- rowSums fine
    Matrix::rowSums(W^2)
  } else {
    # A_new may be sparse (e.g. [I_p | X_grid] for augmented fits).
    # Use Matrix::rowSums to handle both sparse and dense results correctly.
    AnewW <- A_new %*% W              # n_new x q
    Matrix::rowSums(AnewW^2)
  }
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

.rademacher <- function(p, n_probes) {
  matrix(sample(c(-1L, 1L), p * n_probes, replace = TRUE),
         nrow = p, ncol = n_probes)
}

.hutchinson_diag <- function(apply_A, probes) {
  p        <- nrow(probes)
  n_probes <- ncol(probes)
  rowMeans(vapply(seq_len(n_probes), function(i) {
    z <- probes[, i]
    z * apply_A(z)
  }, numeric(p)))
}
