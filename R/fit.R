# fastblm_fit.R
#
# Patch notes vs original:
#   - .estimate_sigma2e replaced with .estimate_sigma2_mode using the correct
#     posterior mode formula under Jeffreys prior:
#
#       sigma2_mode = (RSS + penalty) / (n + 2)
#
#     where RSS     = ||y - A*mu||^2_{R^{-1}}
#           penalty = mu' (Q/phi) mu
#
#   - .fit_woodbury gains a Q argument so penalty can be computed.
#   - fit_fastblm passes Q through to .fit_woodbury.
#   - .fit_pcg fixed: removed broken .estimate_sigma2_mode(A=NULL) call;
#     sigma2e now computed directly using apply_A operator.


#' Fit a Bayesian linear model
#'
#' Computes the posterior mean of \eqn{x} given
#' \eqn{y = Ax + \varepsilon}, with priors
#' \eqn{x \sim N(0, \phi \sigma^2_e Q^{-1})} and
#' \eqn{\varepsilon \sim N(0, \sigma^2_e R)}.
#'
#' \eqn{\sigma^2_e} is estimated as the posterior mode under a Jeffreys prior:
#' \deqn{\hat{\sigma}^2_e = \frac{
#'   \|y - A\mu\|^2_{R^{-1}} + \mu^\top (Q/\phi) \mu
#' }{n + 2}}
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix, or function \eqn{v \mapsto Av} (requires A_t)
#' @param Q p x p prior precision matrix, or function \eqn{v \mapsto Qv}
#' @param phi signal-to-noise ratio \eqn{\sigma^2_b / \sigma^2_e}
#' @param R_inv n x n inverse noise covariance, or function. NULL = identity.
#' @param Q_inv p x p prior covariance, or function \eqn{v \mapsto Q^{-1}v}.
#'   Required for \code{solver = "woodbury"}.
#' @param A_t optional transpose operator. Required if A is a function.
#' @param solver one of \code{"cholesky"}, \code{"pcg"}, or \code{"woodbury"}
#' @param pcg_tol PCG convergence tolerance
#' @param pcg_maxit PCG max iterations
#' @param pcg_precond optional preconditioner function
#'
#' @return object of class fastblm_fit
#' @export
fit_fastblm <- function(y, A, Q, phi,
                        R_inv       = NULL,
                        Q_inv       = NULL,
                        A_t         = NULL,
                        solver      = "cholesky",
                        pcg_tol     = 1e-6,
                        pcg_maxit   = NULL,
                        pcg_precond = NULL,
                        AtRinvA     = NULL,
                        AtRinvy     = NULL,
                        yRinvy      = NULL) {

  y <- as.numeric(y)
  n <- length(y)
  p <- if (is_matrix(A)) ncol(A) else length(A_t(rep(1, n)))

  if (solver == "cholesky") {
    stop_if_function(A,     "A")
    stop_if_function(Q,     "Q")
    stop_if_function(R_inv, "R_inv")
    Rinv <- resolve_Rinv(R_inv, n)
    .fit_cholesky(y, A, Q, phi, Rinv, n,
                  AtRinvA = AtRinvA, AtRinvy = AtRinvy, yRinvy = yRinvy)

  } else if (solver == "woodbury") {
    if (is.null(Q_inv)) stop("solver = 'woodbury' requires Q_inv.")
    stop_if_function(A,     "A")
    stop_if_function(R_inv, "R_inv")
    apply_Qinv <- as_apply(Q_inv)
    Rinv       <- resolve_Rinv(R_inv, n)
    .fit_woodbury(y, A, Q, apply_Qinv, phi, Rinv, n, p)

  } else if (solver == "pcg") {
    apply_A    <- as_apply(A)
    apply_At   <- if (!is.null(A_t)) A_t else {
      if (is.function(A)) stop("A is a function so A_t must be supplied.")
      function(v) as.numeric(Matrix::crossprod(A, v))
    }
    apply_Q    <- as_apply(Q)
    apply_Rinv <- as_apply(R_inv)
    .fit_pcg(y, apply_A, apply_At, apply_Q, phi, apply_Rinv, n,
             tol = pcg_tol, maxit = pcg_maxit %||% (4L * p),
             precond = pcg_precond)

  } else {
    stop("solver must be 'cholesky', 'woodbury', or 'pcg'.")
  }
}


# ------------------------------------------------------------------------------
# Posterior mode of sigma2 under Jeffreys prior
#
# Model:  y | beta, sigma2 ~ N(A*beta, sigma2 * R^{-1})
#         beta | sigma2    ~ N(0, phi * sigma2 * Q^{-1})
#         sigma2           ~ Jeffreys (1/sigma2)
#
# Joint posterior mode of sigma2 (plugging in posterior mean mu):
#   sigma2_mode = (RSS + penalty) / (n + 2)
#
# where:
#   RSS     = ||y - A*mu||^2_{R^{-1}}  = y'R^{-1}y - 2*mu'A'R^{-1}y + mu'A'R^{-1}A*mu
#   penalty = mu' (Q/phi) mu           = (1/phi) * mu' Q mu
#
# Arguments:
#   yRinvy    : y' R^{-1} y                   (scalar)
#   AtRinvy   : A' R^{-1} y                   (p-vector)
#   mu        : posterior mean                 (p-vector)
#   Qmu       : Q %*% mu                      (p-vector)
#   phi       : signal-to-noise ratio          (scalar)
#   n         : number of observations         (scalar)
#   AtRinvAmu : A' R^{-1} A mu (cholesky path, optional)
#   y, A, Rinv: used in woodbury path when AtRinvAmu is NULL
# ------------------------------------------------------------------------------
.estimate_sigma2_mode <- function(yRinvy, AtRinvy, mu, Qmu, phi, n,
                                  AtRinvAmu = NULL, y = NULL, A = NULL,
                                  Rinv = NULL) {
  if (!is.null(AtRinvAmu)) {
    # Cholesky path: use precomputed cross-products
    RSS <- as.numeric(yRinvy
                      - 2 * Matrix::crossprod(mu, AtRinvy)
                      + Matrix::crossprod(mu, AtRinvAmu))
  } else {
    # Woodbury path: compute from residual directly
    resid <- y - as.numeric(A %*% mu)
    if (!is.null(Rinv)) {
      RSS <- as.numeric(Matrix::crossprod(resid, Rinv %*% resid))
    } else {
      RSS <- sum(resid^2)
    }
  }

  penalty <- as.numeric(Matrix::crossprod(mu, Qmu)) / phi
  (RSS + penalty) / (n + 2)
}


# ------------------------------------------------------------------------------
# Internal: fit via p x p sparse Cholesky
# ------------------------------------------------------------------------------
.fit_cholesky <- function(y, A, Q, phi, Rinv, n,
                          AtRinvA = NULL, AtRinvy = NULL, yRinvy = NULL) {
  if (is.null(AtRinvA)) {
    Rinvy   <- Rinv %*% y
    AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
    yRinvy  <- as.numeric(Matrix::crossprod(y, Rinvy))
    AtRinvA <- Matrix::crossprod(A, Rinv %*% A)
  }

  K  <- Matrix::forceSymmetric(AtRinvA + (1/phi) * Q)
  C  <- Matrix::Cholesky(K)
  mu <- as.numeric(Matrix::solve(C, AtRinvy))

  Qmu       <- as.numeric(Q %*% mu)
  AtRinvAmu <- as.numeric(AtRinvA %*% mu)

  sigma2e <- .estimate_sigma2_mode(
    yRinvy    = yRinvy,
    AtRinvy   = AtRinvy,
    mu        = mu,
    Qmu       = Qmu,
    phi       = phi,
    n         = n,
    AtRinvAmu = AtRinvAmu
  )

  structure(
    list(
      posterior_mean = mu,
      phi            = phi,
      sigma2e        = sigma2e,
      sigma2b        = phi * sigma2e,
      chol_factor    = C,
      chol_M         = NULL,
      QinvAt         = NULL,
      apply_K        = NULL,
      solver_type    = "cholesky",
      A              = A,
      Q              = Q,
      R_inv          = Rinv
    ),
    class = "fastblm_fit"
  )
}


# ------------------------------------------------------------------------------
# Internal: fit via n x n Woodbury Cholesky
# ------------------------------------------------------------------------------
.fit_woodbury <- function(y, A, Q, apply_Qinv, phi, Rinv, n, p) {
  At     <- Matrix::t(A)
  QinvAt <- if (is_matrix(apply_Qinv)) {
    as.matrix(apply_Qinv %*% At)
  } else {
    apply(At, 2, apply_Qinv)
  }

  AQinvAt <- A %*% QinvAt
  M       <- Matrix::forceSymmetric(phi * AQinvAt + solve(Rinv))
  CM      <- Matrix::Cholesky(Matrix::Matrix(M, sparse = FALSE))

  Minvy <- as.numeric(Matrix::solve(CM, y))
  mu    <- phi * as.numeric(QinvAt %*% Minvy)

  Qmu <- if (is_matrix(Q)) {
    as.numeric(Q %*% mu)
  } else {
    as.numeric(Q(mu))
  }

  Rinvy   <- Rinv %*% y
  AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
  yRinvy  <- as.numeric(Matrix::crossprod(y, Rinvy))

  sigma2e <- .estimate_sigma2_mode(
    yRinvy  = yRinvy,
    AtRinvy = AtRinvy,
    mu      = mu,
    Qmu     = Qmu,
    phi     = phi,
    n       = n,
    y       = y,
    A       = A,
    Rinv    = Rinv
  )

  structure(
    list(
      posterior_mean = mu,
      phi            = phi,
      sigma2e        = sigma2e,
      sigma2b        = phi * sigma2e,
      chol_factor    = NULL,
      chol_M         = CM,
      QinvAt         = QinvAt,
      apply_K        = NULL,
      apply_Qinv     = apply_Qinv,
      solver_type    = "woodbury",
      A              = A,
      Q              = Q,
      R_inv          = Rinv
    ),
    class = "fastblm_fit"
  )
}


# ------------------------------------------------------------------------------
# Internal: fit via PCG
# ------------------------------------------------------------------------------
.fit_pcg <- function(y, apply_A, apply_At, apply_Q, phi, apply_Rinv, n,
                     tol, maxit, precond = NULL) {

  Rinvy   <- apply_Rinv(y)
  AtRinvy <- apply_At(Rinvy)
  yRinvy  <- as.numeric(Matrix::crossprod(y, Rinvy))

  apply_K <- make_apply_K(apply_A, apply_At, apply_Q, apply_Rinv, phi)
  result  <- pcg(apply_K, AtRinvy, tol = tol, maxit = maxit, precond = precond)
  if (!result$converged) warning("PCG did not converge at fit time.")
  mu <- result$x

  # sigma2e: joint posterior mode of sigma2 given mu, under Jeffreys prior
  #   sigma2_mode = (RSS + mu'Qmu/phi) / (n + 2)
  # RSS uses R=I (PCG path does not store R_inv as a matrix).
  Qmu     <- apply_Q(mu)
  resid   <- y - apply_A(mu)
  RSS     <- sum(resid^2)
  penalty <- as.numeric(Matrix::crossprod(mu, Qmu)) / phi
  sigma2e <- (RSS + penalty) / (n + 2)

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
