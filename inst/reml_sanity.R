# Debug REML formula
# Check what log|Q| looks like vs rho, and compare to a reference implementation
pkgload::load_all(".")
library(Matrix)

make_SAR_Q <- function(p, rho) {
  Q <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                          diagonals = list(rep(-rho, p-1),
                                           c(1, rep(1 + rho^2, p-2), 1),
                                           rep(-rho, p-1)))
  Matrix::forceSymmetric(Q)
}

p   <- 100
n   <- 300
set.seed(42)
A   <- matrix(rnorm(n * p), n, p) / sqrt(p)

rho_true <- 0.7
phi_true <- 8
sigma2e_true <- 0.5

Q_true <- make_SAR_Q(p, rho_true)
x_true <- as.numeric(Matrix::solve(Matrix::Cholesky(Q_true), rnorm(p))) *
  sqrt(phi_true * sigma2e_true)
y      <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

pdf("vignettes/debug_reml.pdf", width = 14, height = 10)
par(mfrow = c(3, 3))

rho_grid <- seq(0.05, 0.95, by = 0.05)

# -----------------------------------------------------------------------
# 1. What does log|Q| look like vs rho?
# -----------------------------------------------------------------------
logdetQ_vals <- sapply(rho_grid, function(rho) {
  Q <- make_SAR_Q(p, rho)
  as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
})

plot(rho_grid, logdetQ_vals, type = "b", pch = 19, col = "#2166ac",
     xlab = "rho", ylab = "log|Q|",
     main = "1. log|Q| vs rho\n(should vary significantly)")
abline(v = rho_true, col = "red", lty = 2)

# -----------------------------------------------------------------------
# 2. What does log|K| look like vs rho at fixed phi?
# -----------------------------------------------------------------------
phi_fixed <- phi_true
logdetK_vals <- sapply(rho_grid, function(rho) {
  Q <- make_SAR_Q(p, rho)
  K <- Matrix::forceSymmetric(Matrix::crossprod(A) + (1/phi_fixed) * Q)
  as.numeric(Matrix::determinant(K, logarithm = TRUE)$modulus)
})

plot(rho_grid, logdetK_vals, type = "b", pch = 19, col = "#e66101",
     xlab = "rho", ylab = "log|K|",
     main = sprintf("2. log|K| vs rho\n(phi fixed = %.1f)", phi_fixed))
abline(v = rho_true, col = "red", lty = 2)

# -----------------------------------------------------------------------
# 3. What does profiled phi look like vs rho?
# -----------------------------------------------------------------------
Rinv    <- resolve_Rinv(NULL, n)
AtRinvy <- as.numeric(Matrix::crossprod(A, y))
yRinvy  <- as.numeric(crossprod(y, y))
probes  <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)

phi_profiled <- sapply(rho_grid, function(rho) {
  Q <- make_SAR_Q(p, rho)
  apply_Q <- as_apply(Q)
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  fastblm:::.profile_phi(
    AtRinvy, yRinvy, A, Rinv, apply_Q, list(Q=Q, Q_matrix=Q),
    p, n, logdetQ, probes, 50L, log(0.01), log(1000),
    "cholesky", 1e-6, 4L*p, rep(0,p)
  )
})

plot(rho_grid, phi_profiled, type = "b", pch = 19, col = "#4dac26",
     xlab = "rho", ylab = "profiled phi",
     main = "3. Profiled phi vs rho")
abline(v = rho_true, col = "red", lty = 2)
abline(h = phi_true, col = "blue", lty = 2)
legend("topright", c(sprintf("True rho=%.1f", rho_true),
                     sprintf("True phi=%.1f", phi_true)),
       col = c("red","blue"), lty = 2, bty = "n")

# -----------------------------------------------------------------------
# 4. Break down each ll component vs rho
# -----------------------------------------------------------------------
ll_components <- t(sapply(seq_along(rho_grid), function(i) {
  rho     <- rho_grid[i]
  Q       <- make_SAR_Q(p, rho)
  apply_Q <- as_apply(Q)
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  phi_i   <- phi_profiled[i]

  res <- fastblm:::.eval_reml_ll(
    phi_i, AtRinvy, yRinvy, A, Rinv, apply_Q, list(Q=Q, Q_matrix=Q),
    p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
  )

  K       <- Matrix::forceSymmetric(Matrix::crossprod(A) + (1/phi_i) * Q)
  logdetK <- as.numeric(Matrix::determinant(K, logarithm = TRUE)$modulus)

  c(ll        = res$ll,
    t1        = -n/2 * log(res$sigma2e),
    t2        = -1/2 * logdetK,
    t3        = -p/2 * log(phi_i),
    t4        = 1/2 * logdetQ,
    sigma2e   = res$sigma2e,
    phi       = phi_i)
}))

plot(rho_grid, ll_components[,"ll"], type = "b", pch = 19, col = "#2166ac",
     xlab = "rho", ylab = "REML ll",
     main = "4. Total REML ll vs rho")
abline(v = rho_true, col = "red", lty = 2)

plot(rho_grid, ll_components[,"t1"], type = "b", pch = 19, col = "#2166ac",
     xlab = "rho", ylab = "value",
     main = "5. -n/2 log(sigma2e) vs rho\n(should be ~flat)")
abline(v = rho_true, col = "red", lty = 2)

plot(rho_grid, ll_components[,"t2"], type = "b", pch = 19, col = "#e66101",
     xlab = "rho", ylab = "value",
     main = "6. -1/2 log|K| vs rho")
abline(v = rho_true, col = "red", lty = 2)

plot(rho_grid, ll_components[,"t3"], type = "b", pch = 19, col = "#4dac26",
     xlab = "rho", ylab = "value",
     main = "7. -p/2 log(phi) vs rho\n(this is the dominant term)")
abline(v = rho_true, col = "red", lty = 2)

plot(rho_grid, ll_components[,"t4"], type = "b", pch = 19, col = "#762a83",
     xlab = "rho", ylab = "value",
     main = "8. +1/2 log|Q| vs rho\n(should provide penalty)")
abline(v = rho_true, col = "red", lty = 2)

# -----------------------------------------------------------------------
# 9. Compare our ll to a reference -- direct marginal ll computation
# The marginal ll of y is:
# log p(y) = -n/2 log(2pi) - 1/2 log|Sigma_y| - 1/2 y' Sigma_y^{-1} y
# where Sigma_y = sigma2e * (I + phi * A Q^{-1} A')
# -----------------------------------------------------------------------
ref_ll <- sapply(seq_along(rho_grid), function(i) {
  rho   <- rho_grid[i]
  Q     <- make_SAR_Q(p, rho)
  phi_i <- phi_profiled[i]

  # sigma2e at this (rho, phi)
  apply_Q <- as_apply(Q)
  res <- fastblm:::.eval_reml_ll(
    phi_i, AtRinvy, yRinvy, A, Rinv, apply_Q, list(Q=Q, Q_matrix=Q),
    p, n, logdetQ_vals[i], probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
  )
  s2 <- res$sigma2e

  # direct marginal ll: y ~ N(0, s2*(I + phi*A Q^{-1} A'))
  Qinv    <- solve(as.matrix(Q))
  Sigma_y <- s2 * (diag(n) + phi_i * A %*% Qinv %*% t(A))
  ld      <- as.numeric(determinant(Sigma_y, logarithm=TRUE)$modulus)
  -n/2 * log(2*pi) - 1/2 * ld - 1/2 * as.numeric(t(y) %*% solve(Sigma_y, y))
})

plot(rho_grid, ref_ll, type = "b", pch = 19, col = "#d01c8b",
     xlab = "rho", ylab = "marginal ll",
     main = "9. Reference marginal ll vs rho\n(ground truth - is this unimodal?)")
abline(v = rho_true, col = "red", lty = 2)

par(mfrow = c(1,1))
dev.off()
cat("Plots saved to vignettes/debug_reml.pdf\n")
