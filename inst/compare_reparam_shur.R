# =============================================================================
# compare_reparam_vs_schur.R
#
# Compares two ways of evaluating the constrained concentrated likelihood
# at fixed (phi, sigma2e):
#
#   1. Reparameterisation: fit reduced model y = A_tilde * gamma + eps
#   2. Schur correction:   unconstrained ll + log|M| correction
#
# Both should give the same answer.
# =============================================================================

set.seed(42)

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

logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)

# Shared
Q_inv <- solve(Q)
K     <- crossprod(A) / sigma2e + (1/(phi*sigma2e)) * Q
K_inv <- solve(K)

# H^{-1} operator (verified correct)
Hinv_apply <- function(v) v - phi^2 * sigma2e * A %*% (K_inv %*% (t(A) %*% v))

# =============================================================================
# Approach 1: Reparameterisation
# =============================================================================
Q_red    <- t(N_mat) %*% Q %*% N_mat
A_tilde  <- A %*% N_mat
K_tilde  <- crossprod(A_tilde) + (1/phi) * Q_red

x_hat_r  <- solve(K_tilde, crossprod(A_tilde, y))
s_r      <- sum(y^2) - as.numeric(crossprod(crossprod(A_tilde, y), x_hat_r))
sigma2e_r <- s_r / n

ll_reparam <- (
  - n/2 * log(sigma2e_r)
  - 1/2 * logdet(K_tilde)
  - p_eff/2 * log(phi)
  + 1/2 * logdet(Q_red)
)

# =============================================================================
# Approach 2: Schur correction
# =============================================================================

# Precompute U (n x r) and Vt (r x n)
CQinvCt <- C %*% Q_inv %*% t(C)                          # r x r
U       <- phi * A %*% Q_inv %*% t(C) %*% solve(CQinvCt) # n x r
Vt      <- C %*% Q_inv %*% t(A)                           # r x n

# M = I_r - Vt H^{-1} U  (r x r, free byproduct)
HinvU   <- Hinv_apply(U)                                   # n x r
M       <- diag(r) - Vt %*% HinvU                         # r x r
Minv    <- solve(M)

# s_c = y' H_c^{-1} y
Hinvy   <- Hinv_apply(y)
s_c     <- as.numeric(crossprod(y,
                                Hinvy + HinvU %*% Minv %*% (Vt %*% Hinvy)))
sigma2e_c <- s_c / n

# Concentrated likelihood:
# l* = -n/2 * log(sigma2e_c) - 1/2*log|K| - p/2*log(phi) + 1/2*log|Q| - 1/2*log|M|
# log|K| and log|Q| computed directly here (would use Lanczos in real code)
ll_schur <- (
  - n/2 * log(sigma2e_c)
  - 1/2 * logdet(K * sigma2e)   # log|K_code| = log|K_H| - p*log(sigma2e)
  # but we want log|K_H| = log|A'A + (1/phi)Q|
  - p/2 * log(phi)
  + 1/2 * logdet(Q)
  - 1/2 * logdet(M)
)

# wait -- need to be careful about which K.
# The concentrated REML ll uses K_H = A'A + (1/phi)Q, not K_code.
# K_code = K_H / sigma2e, so log|K_code| = log|K_H| - p*log(sigma2e)
# => log|K_H| = log|K_code| + p*log(sigma2e)
K_H     <- crossprod(A) + (1/phi) * Q    # = sigma2e * K_code
ll_schur <- (
  - n/2 * log(sigma2e_c)
  - 1/2 * logdet(K_H)
  - p/2 * log(phi)
  + 1/2 * logdet(Q)
  - 1/2 * logdet(M)
)

# =============================================================================
# Results
# =============================================================================
cat("=== Concentrated likelihood comparison ===\n\n")
cat(sprintf("  sigma2e (reparam) : %.8f\n",   sigma2e_r))
cat(sprintf("  sigma2e (schur)   : %.8f\n\n", sigma2e_c))
cat(sprintf("  ll (reparam)      : %.8f\n",   ll_reparam))
cat(sprintf("  ll (schur)        : %.8f\n",   ll_schur))
cat(sprintf("  diff              : %.2e\n\n", ll_schur - ll_reparam))

tol <- 1e-6
if (abs(ll_schur - ll_reparam) < tol) {
  cat("  PASS\n")
} else {
  cat("  FAIL\n")
  cat("\n  Diagnostics:\n")
  cat(sprintf("    log|K_H|        : %.6f\n", logdet(K_H)))
  cat(sprintf("    log|K_tilde|    : %.6f\n", logdet(K_tilde)))
  cat(sprintf("    log|Q|          : %.6f\n", logdet(Q)))
  cat(sprintf("    log|Q_red|      : %.6f\n", logdet(Q_red)))
  cat(sprintf("    log|M|          : %.6f\n", logdet(M)))
  cat(sprintf("    s_r             : %.6f\n", s_r))
  cat(sprintf("    s_c             : %.6f\n", s_c))
}
