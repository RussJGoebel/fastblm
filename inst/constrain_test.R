# constrain() test script
# Verifies correctness of the Schur complement constraint implementation
# Run with: pkgload::load_all("."); source("vignettes/constrain_test.R")

pkgload::load_all(".")
library(Matrix)

cat("=== constrain() correctness tests ===\n\n")

# -----------------------------------------------------------------------
# Setup: small problem so we can compute brute-force ground truth
# -----------------------------------------------------------------------
set.seed(42)
n <- 100
p <- 30
phi <- 5

A <- matrix(rnorm(n * p), n, p)
Q <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                        diagonals = list(rep(-0.5, p-1),
                                         c(1, rep(1.25, p-2), 1),
                                         rep(-0.5, p-1)))
Q <- Matrix::forceSymmetric(Q)

x_true <- rnorm(p)
y      <- as.numeric(A %*% x_true + rnorm(n))

fit <- fit_fastblm(y, A, Q, phi = phi, solver = "cholesky")

# Brute force: form Sigma explicitly
K_mat   <- as.matrix(Matrix::crossprod(A) + (1/phi) * Q)
Sigma   <- solve(K_mat) * fit$sigma2e   # p x p full covariance
mu      <- fit$posterior_mean

cat("Setup: n =", n, "p =", p, "phi =", phi, "\n\n")

# -----------------------------------------------------------------------
# Helper: brute force constrained mean and SE
# -----------------------------------------------------------------------
brute_force_constrain <- function(mu, Sigma, C) {
  SigmaCt    <- Sigma %*% t(C)
  CSigmaCt   <- C %*% SigmaCt
  mu_c       <- mu - SigmaCt %*% solve(CSigmaCt, C %*% mu)
  Sigma_c    <- Sigma - SigmaCt %*% solve(CSigmaCt, t(SigmaCt))
  list(mu = as.numeric(mu_c), se = sqrt(diag(Sigma_c)))
}

# -----------------------------------------------------------------------
# TEST 1: Single sum-to-zero constraint
# -----------------------------------------------------------------------
cat("--- Test 1: single sum-to-zero constraint ---\n")
C1     <- matrix(1, 1, p)
fit_c1 <- constrain(fit, C1)
bf1    <- brute_force_constrain(mu, Sigma, C1)

cat("Constraint satisfied (C mu_c ~0):", round(sum(fit_c1$posterior_mean), 10), "\n")
cat("Max diff mean vs brute force:     ", max(abs(fit_c1$posterior_mean - bf1$mu)), "\n")

se_c1 <- posterior_se(fit_c1)
cat("Max diff SE vs brute force:       ", max(abs(se_c1 - bf1$se)), "\n")
cat("All constrained SE <= unconstrained:", all(se_c1 <= posterior_se(fit) + 1e-10), "\n\n")

# -----------------------------------------------------------------------
# TEST 2: Single contrast constraint x[1] = x[p]
# -----------------------------------------------------------------------
cat("--- Test 2: single contrast x[1] - x[p] = 0 ---\n")
C2     <- matrix(0, 1, p); C2[1] <- 1; C2[p] <- -1
fit_c2 <- constrain(fit, C2)
bf2    <- brute_force_constrain(mu, Sigma, C2)

cat("Constraint satisfied (x[1]-x[p]):", round(fit_c2$posterior_mean[1] - fit_c2$posterior_mean[p], 10), "\n")
cat("Max diff mean vs brute force:    ", max(abs(fit_c2$posterior_mean - bf2$mu)), "\n")

se_c2 <- posterior_se(fit_c2)
cat("Max diff SE vs brute force:      ", max(abs(se_c2 - bf2$se)), "\n\n")

# -----------------------------------------------------------------------
# TEST 3: Multiple constraints
# -----------------------------------------------------------------------
cat("--- Test 3: two constraints simultaneously ---\n")
C3     <- rbind(C1, C2)
fit_c3 <- constrain(fit, C3)
bf3    <- brute_force_constrain(mu, Sigma, C3)

residuals <- as.numeric(C3 %*% fit_c3$posterior_mean)
cat("Constraint residuals:", round(residuals, 10), "\n")
cat("Max diff mean vs brute force:", max(abs(fit_c3$posterior_mean - bf3$mu)), "\n")

se_c3 <- posterior_se(fit_c3)
cat("Max diff SE vs brute force:  ", max(abs(se_c3 - bf3$se)), "\n\n")

# -----------------------------------------------------------------------
# TEST 4: Idempotency -- applying same constraint twice = once
# -----------------------------------------------------------------------
cat("--- Test 4: idempotency ---\n")
fit_c1_twice <- constrain(fit_c1, C1)
cat("Max diff mean (once vs twice):", max(abs(fit_c1$posterior_mean - fit_c1_twice$posterior_mean)), "\n")
se_once  <- posterior_se(fit_c1)
se_twice <- posterior_se(fit_c1_twice)
cat("Max diff SE (once vs twice):  ", max(abs(se_once - se_twice)), "\n\n")

# -----------------------------------------------------------------------
# TEST 5: Correction lies in row space of C
# -----------------------------------------------------------------------
cat("--- Test 5: correction in row space of C ---\n")
correction <- fit_c1$posterior_mean - mu
# project onto null space of C -- should be zero
null_proj <- correction - t(C1) %*% solve(C1 %*% t(C1), C1 %*% correction)
cat("Max component orthogonal to row(C):", max(abs(null_proj)), "\n\n")

# -----------------------------------------------------------------------
# TEST 6: consistency across solvers
# -----------------------------------------------------------------------
cat("--- Test 6: Cholesky vs PCG constrained mean ---\n")
fit_pcg <- fit_fastblm(y, A, Q, phi = phi, solver = "pcg")
fit_pcg_c <- constrain(fit_pcg, C1)
cat("Max diff constrained mean (Cholesky vs PCG):",
    max(abs(fit_c1$posterior_mean - fit_pcg_c$posterior_mean)), "\n\n")

# -----------------------------------------------------------------------
# TEST 7: Woodbury constrained mean
# -----------------------------------------------------------------------
cat("--- Test 7: Woodbury constrained mean ---\n")
CQ     <- Matrix::Cholesky(Q)
apply_Qinv <- function(v) as.numeric(Matrix::solve(CQ, v))
fit_wb <- fit_fastblm(y, A, Q, phi = phi, Q_inv = apply_Qinv, solver = "woodbury")
fit_wb_c <- constrain(fit_wb, C1)
cat("Max diff constrained mean (Cholesky vs Woodbury):",
    max(abs(fit_c1$posterior_mean - fit_wb_c$posterior_mean)), "\n\n")

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------
pdf("vignettes/constrain_test.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

# Mean comparison: unconstrained vs constrained
plot(mu, fit_c1$posterior_mean,
     xlab = "Unconstrained mean", ylab = "Constrained mean",
     main = "Test 1: Sum-to-zero constraint\nposterior mean",
     pch = 19, col = "#2166ac", cex = 0.8)
abline(0, 1, col = "red", lty = 2)

# SE comparison
se_unc <- posterior_se(fit)
plot(se_unc, se_c1,
     xlab = "Unconstrained SE", ylab = "Constrained SE",
     main = "Test 1: Sum-to-zero constraint\nposterior SE",
     pch = 19, col = "#2166ac", cex = 0.8)
abline(0, 1, col = "red", lty = 2)

# Brute force vs package SE
plot(bf1$se, se_c1,
     xlab = "Brute force SE", ylab = "Package SE",
     main = "Test 1: SE brute force vs package",
     pch = 19, col = "#4dac26", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", sprintf("max diff = %.2e", max(abs(se_c1 - bf1$se))), bty = "n")

# Two constraints SE
plot(se_unc, se_c3,
     xlab = "Unconstrained SE", ylab = "2-constraint SE",
     main = "Test 3: Two constraints\nposterior SE",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)

# Brute force vs package, two constraints
plot(bf3$se, se_c3,
     xlab = "Brute force SE", ylab = "Package SE",
     main = "Test 3: Two constraints\nSE brute force vs package",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", sprintf("max diff = %.2e", max(abs(se_c3 - bf3$se))), bty = "n")

# Solver consistency
plot(fit_c1$posterior_mean, fit_pcg_c$posterior_mean,
     xlab = "Cholesky constrained mean", ylab = "PCG constrained mean",
     main = "Test 6: Solver consistency",
     pch = 19, col = "#e66101", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", sprintf("max diff = %.2e",
                          max(abs(fit_c1$posterior_mean - fit_pcg_c$posterior_mean))), bty = "n")

par(mfrow = c(1, 1))
dev.off()

cat("=== All constrain tests complete ===\n")
cat("Plots saved to vignettes/constrain_test.pdf\n")
