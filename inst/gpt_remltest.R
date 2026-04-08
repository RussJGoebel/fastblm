# Quick sanity test for exact REML against a dense reference
#
# This script compares the rewritten REML code to a brute-force dense
# calculation on a tiny toy problem. The rewritten code drops the constant
# term -1/2 log|R| from the log-likelihood when R is fixed, so the script
# checks equality up to that constant.

library(Matrix)

# -----------------------------------------------------------------------------
# Dense brute-force REML log-likelihood
# -----------------------------------------------------------------------------
#
# Model:
#   y = A x + X beta + e
#   e ~ N(0, sigma2e * R)
#   x ~ N(0, sigma2e * phi * Q^{-1})
#
# Optional exact constraint:
#   C x = 0
#
# Marginal covariance:
#   H  = R + phi A Q^{-1} A'
# or, with constraint,
#   Hc = R + phi A S A'
# where
#   S = Q^{-1} - Q^{-1} C' (C Q^{-1} C')^{-1} C Q^{-1}
#
# Restricted log-likelihood up to constants:
#   -1/2 log|Hc| - 1/2 log|X'Hc^{-1}X| - (n_eff/2) log(quad / n_eff)
#
dense_reml_ll <- function(y, A, Q, phi, X_fixed = NULL, R = NULL, C = NULL) {
  y <- as.numeric(y)
  A <- as.matrix(A)
  Q <- as.matrix(Q)

  n <- length(y)

  if (is.null(R)) {
    R <- diag(n)
  } else {
    R <- as.matrix(R)
  }

  Q_inv <- solve(Q)

  if (is.null(C)) {
    S <- Q_inv
  } else {
    C <- as.matrix(C)
    CQinvCt <- C %*% Q_inv %*% t(C)
    S <- Q_inv - Q_inv %*% t(C) %*% solve(CQinvCt) %*% C %*% Q_inv
  }

  H <- R + phi * A %*% S %*% t(A)
  H <- 0.5 * (H + t(H))
  H_inv <- solve(H)

  if (is.null(X_fixed)) {
    n_eff <- n
    logdet_XHX <- 0
    quad <- drop(crossprod(y, H_inv %*% y))
  } else {
    X_fixed <- as.matrix(X_fixed)

    XtHinvX <- crossprod(X_fixed, H_inv %*% X_fixed)
    XtHinvX <- 0.5 * (XtHinvX + t(XtHinvX))

    b <- as.numeric(crossprod(X_fixed, H_inv %*% y))
    alpha <- solve(XtHinvX, b)

    quad <- drop(crossprod(y, H_inv %*% y)) - drop(crossprod(b, alpha))
    n_eff <- n - qr(X_fixed)$rank
    logdet_XHX <- as.numeric(determinant(XtHinvX, logarithm = TRUE)$modulus)
  }

  logdet_H <- as.numeric(determinant(H, logarithm = TRUE)$modulus)
  logdet_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)

  list(
    ll_full = -0.5 * logdet_H -
      0.5 * logdet_XHX -
      (n_eff / 2) * log(quad / n_eff),
    ll_drop_R = -0.5 * (logdet_H - logdet_R) -
      0.5 * logdet_XHX -
      (n_eff / 2) * log(quad / n_eff),
    quad = quad,
    n_eff = n_eff,
    logdet_R = logdet_R
  )
}

# -----------------------------------------------------------------------------
# Tiny test problem
# -----------------------------------------------------------------------------
set.seed(123)

n <- 6
p <- 4

A <- matrix(rnorm(n * p), n, p)

# Positive definite Q
B <- matrix(rnorm(p * p), p, p)
Q <- crossprod(B) + diag(0.5, p)

# Fixed effects
X_fixed <- cbind(1, seq_len(n))

# Observation covariance
R <- diag(c(1.0, 1.2, 0.8, 1.1, 0.9, 1.3))
R_inv <- solve(R)

# Constraint: first two latent coefficients sum to zero
C <- matrix(c(1, 1, 0, 0), nrow = 1)

# Simulate a response
phi_true <- 2.0
sigma2e_true <- 1.5

x_true <- drop(chol(sigma2e_true * phi_true * solve(Q)) %*% rnorm(p))
beta_true <- c(0.3, -0.2)
e_true <- drop(chol(sigma2e_true * R) %*% rnorm(n))
y <- drop(A %*% x_true + X_fixed %*% beta_true + e_true)

# -----------------------------------------------------------------------------
# Prior object for tune_reml() internals
# -----------------------------------------------------------------------------
Q_fun_simple <- function(theta) {
  list(
    Q = Q,
    Q_matrix = Q,
    log_det_Q = as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  )
}

# -----------------------------------------------------------------------------
# Wrapper to evaluate the rewritten REML objective at fixed phi
# -----------------------------------------------------------------------------
eval_new_ll <- function(y, A, phi, X_fixed = NULL, R_inv = NULL, C = NULL) {
  n <- length(y)
  p <- ncol(A)

  Rinv_obj <- .make_Rinv(R_inv, n)
  apply_Rinv <- Rinv_obj$apply

  if (!is.null(Rinv_obj$matrix)) {
    AtRinvA_matrix <- crossprod(A, Rinv_obj$matrix %*% A)
  } else {
    RinvA <- apply(A, 2, apply_Rinv)
    RinvA <- as.matrix(RinvA)
    AtRinvA_matrix <- crossprod(A, RinvA)
  }

  AtRinvA_apply <- function(v) as.numeric(AtRinvA_matrix %*% v)

  prior <- Q_fun_simple(numeric(0))
  apply_Q <- .as_apply(prior$Q)

  set.seed(1)
  probes <- matrix(sample(c(-1L, 1L), p * 10, replace = TRUE), nrow = p)

  .eval_reml_ll(
    y = y,
    A = A,
    X_fixed = X_fixed,
    phi = phi,
    prior = prior,
    apply_Q = apply_Q,
    logdet_Q = prior$log_det_Q,
    apply_Rinv = apply_Rinv,
    AtRinvA_apply = AtRinvA_apply,
    AtRinvA_matrix = AtRinvA_matrix,
    probes = probes,
    n_lanczos_steps = 10L,
    logdet_method = "cholesky",
    pcg_tol = 1e-10,
    pcg_maxit = 1000L,
    constraint_matrix = C
  )
}

# -----------------------------------------------------------------------------
# Run four comparisons
# -----------------------------------------------------------------------------
cases <- list(
  list(name = "no_fixed_no_constraint", X = NULL,    C = NULL),
  list(name = "fixed_no_constraint",    X = X_fixed, C = NULL),
  list(name = "no_fixed_constraint",    X = NULL,    C = C),
  list(name = "fixed_constraint",       X = X_fixed, C = C)
)

results <- lapply(cases, function(case) {
  dense <- dense_reml_ll(
    y = y,
    A = A,
    Q = Q,
    phi = phi_true,
    X_fixed = case$X,
    R = R,
    C = case$C
  )

  new <- eval_new_ll(
    y = y,
    A = A,
    phi = phi_true,
    X_fixed = case$X,
    R_inv = R_inv,
    C = case$C
  )

  data.frame(
    case = case$name,
    dense_ll_full = dense$ll_full,
    dense_ll_drop_R = dense$ll_drop_R,
    new_ll = new$ll,
    full_diff = abs(dense$ll_full - new$ll),
    drop_R_diff = abs(dense$ll_drop_R - new$ll),
    dense_sigma2e = dense$quad / dense$n_eff,
    new_sigma2e = new$sigma2e,
    sigma2e_diff = abs(dense$quad / dense$n_eff - new$sigma2e)
  )
})

results_df <- do.call(rbind, results)
print(results_df)

# -----------------------------------------------------------------------------
# Check expected offset
# -----------------------------------------------------------------------------
logdet_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
expected_offset <- 0.5 * logdet_R

cat("\nExpected dropped constant 0.5 * log|R|:\n")
print(expected_offset)

cat("\nObserved full-difference values:\n")
print(results_df$full_diff)

cat("\nObserved difference after removing log|R| constant:\n")
print(results_df$drop_R_diff)

# -----------------------------------------------------------------------------
# Hard checks
# -----------------------------------------------------------------------------
tol <- 1e-5

if (any(abs(results_df$full_diff - expected_offset) > tol)) {
  stop("Full log-likelihood mismatch is not just the expected 0.5 * log|R| constant.")
}

if (any(results_df$drop_R_diff > tol)) {
  stop("Log-likelihood mismatch exceeds tolerance after removing the log|R| constant.")
}

if (any(results_df$sigma2e_diff > tol)) {
  stop("Profiled sigma2e mismatch exceeds tolerance.")
}

cat("\nAll REML sanity checks passed.\n")
