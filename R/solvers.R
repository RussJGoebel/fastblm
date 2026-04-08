#' PCG solver
#'
#' Solves \eqn{Ax = b} iteratively via preconditioned conjugate gradient.
#'
#' @param apply_A function \eqn{v \mapsto Av}
#' @param b right hand side vector
#' @param x0 initial guess, default zero
#' @param precond optional preconditioner function \eqn{v \mapsto M^{-1} v}
#' @param tol convergence tolerance
#' @param maxit maximum iterations
#'
#' @return list with elements \code{x} (solution), \code{converged} (logical),
#'   and \code{iter} (number of iterations taken)
#' @export
pcg <- function(apply_A, b, x0 = NULL, precond = NULL, tol = 1e-6, maxit = NULL) {
  p     <- length(b)
  maxit <- maxit %||% (4L * p)
  x     <- x0 %||% rep(0, p)

  precond <- precond %||% function(v) v

  r  <- b - apply_A(x)
  z  <- precond(r)
  d  <- z
  rz <- as.numeric(crossprod(r, z))

  for (i in seq_len(maxit)) {
    Ad    <- apply_A(d)
    dAd   <- as.numeric(crossprod(d, Ad))
    if (dAd <= 0) stop("PCG: matrix appears non-positive-definite.")
    alpha <- rz / dAd
    x     <- x + alpha * d
    r     <- r - alpha * Ad
    if (sqrt(sum(r^2)) < tol) return(list(x = x, converged = TRUE, iter = i))
    z     <- precond(r)
    rz_new <- as.numeric(crossprod(r, z))
    beta  <- rz_new / rz
    d     <- z + beta * d
    rz    <- rz_new
  }

  warning("PCG did not converge in ", maxit, " iterations.")
  list(x = x, converged = FALSE, iter = maxit)
}

#' Cholesky log determinant
#'
#' Computes \eqn{\log |A|} from a Cholesky factor.
#'
#' @param C CHMfactor from \code{Matrix::Cholesky}
#'
#' @return numeric scalar
#' @export
chol_logdet <- function(C) {
  2 * sum(log(Matrix::diag(Matrix::expand(C)$L)))
}

#' Stochastic log determinant via Lanczos quadrature
#'
#' Estimates \eqn{\log |A|} using stochastic Lanczos quadrature with
#' Rademacher probe vectors.
#'
#' @param apply_A function \eqn{v \mapsto Av} (must be symmetric positive definite)
#' @param probes matrix of Rademacher probe vectors (p x n_probes)
#' @param n_steps number of Lanczos steps
#'
#' @return numeric scalar estimate of \eqn{\log |A|}
#' @export
lanczos_logdet <- function(apply_A, probes, n_steps = 50L) {
  p        <- nrow(probes)
  n_probes <- ncol(probes)

  estimates <- vapply(seq_len(n_probes), function(i) {
    lanczos_quadrature(apply_A, probes[, i], log, n_steps)$estimate
  }, numeric(1L))

  mean(estimates)
}

#' Single Lanczos quadrature estimate for f(A) applied to probe vector
#'
#' Uses Lanczos tridiagonalization and Gauss quadrature to estimate
#' \eqn{v^\top f(A) v} for a scalar function \eqn{f}.
#'
#' @param apply_A function \eqn{v \mapsto Av}
#' @param probe probe vector
#' @param f scalar function to apply to eigenvalues
#' @param n_steps number of Lanczos steps
#'
#' @return list with element \code{estimate}
#' @export
lanczos_quadrature <- function(apply_A, probe, f, n_steps) {
  p      <- length(probe)
  n_steps <- min(n_steps, p)

  # Lanczos tridiagonalization
  Q      <- matrix(0, p, n_steps + 1)
  alpha  <- numeric(n_steps)
  beta   <- numeric(n_steps)

  Q[, 1] <- probe / sqrt(sum(probe^2))

  for (j in seq_len(n_steps)) {
    z         <- apply_A(Q[, j])
    alpha[j]  <- as.numeric(crossprod(Q[, j], z))
    z         <- z - alpha[j] * Q[, j]
    if (j > 1) z <- z - beta[j-1] * Q[, j-1]
    beta[j]   <- sqrt(sum(z^2))
    if (beta[j] < 1e-12) { n_steps <- j; break }
    Q[, j+1]  <- z / beta[j]
  }

  # Tridiagonal eigendecomposition
  alpha <- alpha[seq_len(n_steps)]
  beta  <- beta[seq_len(n_steps - 1)]
  T_mat <- diag(alpha) +
    if (n_steps > 1) {
      diag(beta, n_steps, n_steps - 1) %*% diag(1, n_steps - 1, n_steps) +
        diag(1, n_steps, n_steps - 1) %*% diag(beta, n_steps - 1, n_steps)
    } else matrix(0, 1, 1)

  eig    <- eigen(T_mat, symmetric = TRUE)
  tau    <- eig$values
  e1     <- eig$vectors[1, ]
  norm2  <- sum(probe^2)

  estimate <- norm2 * sum(e1^2 * f(pmax(tau, .Machine$double.eps)))
  list(estimate = estimate)
}
