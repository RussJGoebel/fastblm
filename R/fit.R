#' Fit a Bayesian linear model
#'
#' Computes the posterior mean of x given y = Ax + epsilon,
#' x ~ N(0, phi * sigma2e * Q^{-1}), epsilon ~ N(0, sigma2e * R).
#'
#' sigma2e is estimated as the posterior mode given phi:
#' sigma2e_hat = (y' R^{-1} y - mu' A' R^{-1} y) / n
#'
#' Three solvers are available:
#'   "cholesky" -- p x p Cholesky factor. Best when p << n.
#'   "pcg"      -- iterative PCG. Best when p >> n and Q^{-1} is not cheap.
#'   "woodbury" -- n x n Cholesky via Woodbury. Best when p >> n and Q^{-1} is cheap.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix, or function v -> Av (requires A_t)
#' @param Q p x p prior precision matrix, or function v -> Qv
#' @param phi signal-to-noise ratio sigma2b / sigma2e
#' @param R_inv n x n inverse noise covariance, or function v -> R_inv v. NULL = identity.
#' @param Q_inv p x p prior covariance matrix, or function v -> Q^{-1}v. Required for solver = "woodbury".
#' @param A_t optional transpose operator function v -> A'v. Required if A is a function.
#' @param solver "cholesky", "pcg", or "woodbury"
#' @param pcg_tol PCG convergence tolerance
#' @param pcg_maxit PCG max iterations
#'
#' @return object of class fastblm_fit
#' @export
fit_fastblm <- function(y, A, Q, phi,
                        R_inv     = NULL,
                        Q_inv     = NULL,
                        A_t       = NULL,
                        solver    = "cholesky",
                        pcg_tol   = 1e-6,
                        pcg_maxit = NULL) {

  y <- as.numeric(y)
  n <- length(y)
  p <- if (is_matrix(A)) ncol(A) else length(A_t(rep(1, n)))

  if (solver == "cholesky") {
    stop_if_function(A,     "A")
    stop_if_function(Q,     "Q")
    stop_if_function(R_inv, "R_inv")
    Rinv <- resolve_Rinv(R_inv, n)
    .fit_cholesky(y, A, Q, phi, Rinv, n)

  } else if (solver == "woodbury") {
    if (is.null(Q_inv)) stop("solver = 'woodbury' requires Q_inv.")
    stop_if_function(A,     "A")
    stop_if_function(R_inv, "R_inv")
    apply_Qinv <- as_apply(Q_inv)
    Rinv       <- resolve_Rinv(R_inv, n)
    .fit_woodbury(y, A, apply_Qinv, phi, Rinv, n, p)

  } else if (solver == "pcg") {
    apply_A    <- as_apply(A)
    apply_At   <- if (!is.null(A_t)) A_t else {
      if (is.function(A)) stop("A is a function so A_t must be supplied.")
      function(v) as.numeric(Matrix::crossprod(A, v))
    }
    apply_Q    <- as_apply(Q)
    apply_Rinv <- as_apply(R_inv)
    .fit_pcg(y, apply_A, apply_At, apply_Q, phi, apply_Rinv, n,
             tol = pcg_tol, maxit = pcg_maxit %||% (4L * p))

  } else {
    stop("solver must be 'cholesky', 'woodbury', or 'pcg'.")
  }
}

# Estimate sigma2e as posterior mode given phi and posterior mean
# sigma2e = (y' R^{-1} y - mu' A' R^{-1} y) / n
.estimate_sigma2e <- function(yRinvy, AtRinvy, mu, n) {
  as.numeric((yRinvy - crossprod(mu, AtRinvy)) / n)
}

# Internal: fit via p x p sparse Cholesky
.fit_cholesky <- function(y, A, Q, phi, Rinv, n) {
  Rinvy   <- Rinv %*% y
  AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
  yRinvy  <- as.numeric(crossprod(y, Rinvy))

  AtRinvA <- Matrix::crossprod(A, Rinv %*% A)
  K       <- Matrix::forceSymmetric(AtRinvA + (1/phi) * Q)
  C       <- Matrix::Cholesky(K)
  mu      <- as.numeric(Matrix::solve(C, AtRinvy))
  sigma2e <- .estimate_sigma2e(yRinvy, AtRinvy, mu, n)

  structure(
    list(
      posterior_mean = mu,
      phi            = phi,
      sigma2e        = sigma2e,
      sigma2b        = phi * sigma2e,
      chol_factor    = C,          # p x p CHMfactor
      chol_M         = NULL,       # n x n CHMfactor (woodbury only)
      QinvAt         = NULL,       # p x n matrix (woodbury only)
      apply_K        = NULL,
      solver_type    = "cholesky",
      A              = A,
      Q              = Q,
      R_inv          = Rinv
    ),
    class = "fastblm_fit"
  )
}

# Internal: fit via n x n Woodbury Cholesky
# Uses: x = phi Q^{-1} A' M^{-1} y  where M = phi A Q^{-1} A' + R
# Cheap when p >> n and Q^{-1} is easy to apply
.fit_woodbury <- function(y, A, apply_Qinv, phi, Rinv, n, p) {
  # form Q^{-1} A' -- p x n matrix, requires n applications of Q^{-1}
  At     <- Matrix::t(A)
  QinvAt <- apply(At, 2, apply_Qinv)   # p x n

  # M = phi A Q^{-1} A' + R  -- n x n
  AQinvAt <- A %*% QinvAt              # n x n
  M       <- Matrix::forceSymmetric(phi * AQinvAt + solve(Rinv))
  CM      <- Matrix::Cholesky(Matrix::Matrix(M, sparse = FALSE))

  # posterior mean: x = phi Q^{-1} A' M^{-1} y
  Minvy   <- as.numeric(Matrix::solve(CM, y))
  mu      <- phi * as.numeric(QinvAt %*% Minvy)

  # sigma2e
  Rinvy   <- Rinv %*% y
  AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
  yRinvy  <- as.numeric(crossprod(y, Rinvy))
  sigma2e <- .estimate_sigma2e(yRinvy, AtRinvy, mu, n)

  structure(
    list(
      posterior_mean = mu,
      phi            = phi,
      sigma2e        = sigma2e,
      sigma2b        = phi * sigma2e,
      chol_factor    = NULL,
      chol_M         = CM,         # n x n CHMfactor — reused in posterior_se
      QinvAt         = QinvAt,     # p x n — reused in posterior_se
      apply_K        = NULL,
      apply_Qinv     = apply_Qinv, # reused for diag(Q^{-1}) in posterior_se
      solver_type    = "woodbury",
      A              = A,
      Q              = NULL,
      R_inv          = Rinv
    ),
    class = "fastblm_fit"
  )
}

# Internal: fit via PCG
.fit_pcg <- function(y, apply_A, apply_At, apply_Q, phi, apply_Rinv, n,
                     tol, maxit) {
  Rinvy   <- apply_Rinv(y)
  AtRinvy <- apply_At(Rinvy)
  yRinvy  <- as.numeric(crossprod(y, Rinvy))

  apply_K <- make_apply_K(apply_A, apply_At, apply_Q, apply_Rinv, phi)
  result  <- pcg(apply_K, AtRinvy, tol = tol, maxit = maxit)
  if (!result$converged) warning("PCG did not converge at fit time.")
  mu      <- result$x
  sigma2e <- .estimate_sigma2e(yRinvy, AtRinvy, mu, n)

  structure(
    list(
      posterior_mean = mu,
      phi            = phi,
      sigma2e        = sigma2e,
      sigma2b        = phi * sigma2e,
      chol_factor    = NULL,
      chol_M         = NULL,
      QinvAt         = NULL,
      apply_K        = apply_K,
      solver_type    = "pcg",
      A              = NULL,
      Q              = NULL,
      R_inv          = NULL
    ),
    class = "fastblm_fit"
  )
}

#' @export
print.fastblm_fit <- function(x, ...) {
  cat("fastblm_fit\n")
  cat(sprintf("  p        : %d coefficients\n", length(x$posterior_mean)))
  cat(sprintf("  solver   : %s\n", x$solver_type))
  cat(sprintf("  phi      : %.4g\n", x$phi))
  cat(sprintf("  sigma2e  : %.4g\n", x$sigma2e))
  cat(sprintf("  sigma2b  : %.4g\n", x$sigma2b))
  invisible(x)
}
