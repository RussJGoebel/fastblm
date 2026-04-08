# =============================================================================
# verify_constrained_ll.R  (v2)
#
# Evaluates log p(y | Cx=0) at FIXED (phi, sigma2e) -- no profiling.
# All three approaches should give the same number.
#
#   1. Monte Carlo: sample gamma from prior, average p(y | N*gamma)
#   2. Reparameterisation: analytic Gaussian integral in reduced space
#   3. Correction: unconstrained analytic integral + correction terms
# =============================================================================

set.seed(42)
library(Matrix)

n     <- 20
p     <- 6
r     <- 2
p_eff <- p - r

A <- matrix(rnorm(n * p), n, p)

Q <- diag(p) + 0.3 * crossprod(matrix(rnorm(p * p), p, p)) / p
Q <- (Q + t(Q)) / 2
Q <- Q + p * diag(p)

phi     <- 2.0
sigma2e <- 0.5

C <- matrix(rnorm(r * p), r, p)
while (qr(C)$rank < r) C <- matrix(rnorm(r * p), r, p)

qr_C  <- qr(t(C))
N_mat <- qr.Q(qr_C, complete = TRUE)[, (r + 1):p, drop = FALSE]

x_true <- N_mat %*% rnorm(p_eff)
y      <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e)))

cat(sprintf("Setup: n=%d, p=%d, r=%d, p_eff=%d, phi=%.2f, sigma2e=%.2f\n\n",
            n, p, r, p_eff, phi, sigma2e))

logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)

Q_inv <- solve(Q)
K     <- crossprod(A) / sigma2e + (1/(phi*sigma2e)) * Q
K_inv <- solve(K)

# --- Approach 1: Monte Carlo ---
cat("Approach 1: Monte Carlo\n")
Q_red       <- t(N_mat) %*% Q %*% N_mat
Sigma_gamma <- phi * sigma2e * solve(Q_red)
L_mc        <- chol(Sigma_gamma)
n_mc        <- 1000000L
Gamma       <- L_mc %*% matrix(rnorm(p_eff * n_mc), p_eff, n_mc)
resid       <- y - A %*% (N_mat %*% Gamma)
log_lik     <- -n/2 * log(2*pi*sigma2e) - colSums(resid^2) / (2*sigma2e)
max_ll      <- max(log_lik)
ll_mc       <- max_ll + log(mean(exp(log_lik - max_ll)))
cat(sprintf("  log p(y | Cx=0) = %.6f\n\n", ll_mc))

# --- Approach 2: Reparameterisation (analytic) ---
cat("Approach 2: Reparameterisation (analytic)\n")
A_tilde    <- A %*% N_mat
H_tilde    <- diag(n) + phi * A_tilde %*% solve(Q_red) %*% t(A_tilde)
ll_reparam <- (
  - n/2 * log(2*pi)
  - n/2 * log(sigma2e)
  - 1/2 * logdet(H_tilde)
  - 1/(2*sigma2e) * as.numeric(crossprod(y, solve(H_tilde, y)))
)
cat(sprintf("  log p(y | Cx=0) = %.6f\n\n", ll_reparam))

# --- Approach 3: Bayes correction ---
cat("Approach 3: Correction (Bayes)\n")
H       <- diag(n) + phi * A %*% Q_inv %*% t(A)
H_inv   <- solve(H)
ll_unc  <- (
  - n/2 * log(2*pi)
  - n/2 * log(sigma2e)
  - 1/2 * logdet(H)
  - 1/(2*sigma2e) * as.numeric(crossprod(y, H_inv %*% y))
)
x_hat       <- phi * Q_inv %*% t(A) %*% H_inv %*% y
Cx_hat      <- as.numeric(C %*% x_hat)
CKinvCt     <- sigma2e * C %*% K_inv %*% t(C)
CQinvCt     <- phi * sigma2e * C %*% Q_inv %*% t(C)
log_post_at_0 <- (
  - r/2 * log(2*pi)
  - 1/2 * logdet(CKinvCt)
  - 1/2 * as.numeric(crossprod(Cx_hat, solve(CKinvCt, Cx_hat)))
)
log_prior_at_0 <- (
  - r/2 * log(2*pi)
  - 1/2 * logdet(CQinvCt)
)
ll_correction <- ll_unc + log_post_at_0 - log_prior_at_0
cat(sprintf("  ll_unconstrained  = %.6f\n", ll_unc))
cat(sprintf("  log p(Cx=0|y)     = %.6f\n", log_post_at_0))
cat(sprintf("  log p(Cx=0)       = %.6f\n", log_prior_at_0))
cat(sprintf("  log p(y | Cx=0)   = %.6f\n\n", ll_correction))

cat("=== Summary ===\n")
cat(sprintf("  Monte Carlo:           %.6f\n", ll_mc))
cat(sprintf("  Reparameterisation:    %.6f\n", ll_reparam))
cat(sprintf("  Correction:            %.6f\n", ll_correction))
cat(sprintf("\n  Reparam   vs MC:       diff = %.6f\n", ll_reparam - ll_mc))
cat(sprintf("  Correction vs Reparam: diff = %.6f\n", ll_correction - ll_reparam))
cat(sprintf("  Correction vs MC:      diff = %.6f\n", ll_correction - ll_mc))

if (abs(ll_correction - ll_reparam) < 0.05) {
  cat("\n  PASS: correction agrees with reparameterisation\n")
} else {
  cat("\n  FAIL: approaches disagree\n")
}
