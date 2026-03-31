# fastblm sanity checks
# Run with: source("vignettes/sanity_checks.R")
# or: pkgload::load_all("."); source("vignettes/sanity_checks.R")

library(Matrix)
pkgload::load_all(".")

pdf("vignettes/sanity_checks.pdf", width = 10, height = 8)

cat("=== fastblm sanity checks ===\n\n")

# -----------------------------------------------------------------------
# Simulate data
# y = A x + epsilon
# x ~ N(0, phi * sigma2e * Q^{-1})
# epsilon ~ N(0, sigma2e * I)
# -----------------------------------------------------------------------
set.seed(42)
n   <- 200
p   <- 50
phi     <- 5
sigma2e <- 1
sigma2b <- phi * sigma2e

# A: simple random design matrix
A <- matrix(rnorm(n * p), n, p)

# Q: diagonal precision (iid prior on x)
Q <- Matrix::Diagonal(p, 1)

# true x
x_true <- rnorm(p, sd = sqrt(sigma2b))

# response
y <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e)))

cat("Simulated data: n =", n, "p =", p, "\n")
cat("True phi =", phi, "  sigma2e =", sigma2e, "\n\n")

# -----------------------------------------------------------------------
# CHECK 1: Cholesky fit, R_inv = identity (default)
# -----------------------------------------------------------------------
cat("--- Check 1: Cholesky fit, R_inv = NULL (identity) ---\n")
fit_chol <- fit_fastblm(y, A, Q, phi = phi, solver = "cholesky")
print(fit_chol)
cat("posterior_mean range:", range(fit_chol$posterior_mean), "\n")
cat("correlation with x_true:", cor(fit_chol$posterior_mean, x_true), "\n\n")

# Plot: posterior mean vs true x
plot(x_true, fit_chol$posterior_mean,
     xlab = "True x", ylab = "Posterior mean",
     main = "Check 1: Cholesky fit — posterior mean vs true x",
     pch = 19, col = "#2166ac", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", legend = sprintf("cor = %.3f", cor(fit_chol$posterior_mean, x_true)),
       bty = "n")

# -----------------------------------------------------------------------
# CHECK 2: PCG fit, should give same posterior mean as Cholesky
# -----------------------------------------------------------------------
cat("--- Check 2: PCG fit ---\n")
fit_pcg <- fit_fastblm(y, A, Q, phi = phi, solver = "pcg")
print(fit_pcg)
cat("Max diff vs Cholesky:", max(abs(fit_pcg$posterior_mean - fit_chol$posterior_mean)), "\n\n")

# Plot: PCG vs Cholesky posterior means
plot(fit_chol$posterior_mean, fit_pcg$posterior_mean,
     xlab = "Cholesky posterior mean", ylab = "PCG posterior mean",
     main = "Check 2: PCG vs Cholesky posterior means",
     pch = 19, col = "#4dac26", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", legend = sprintf("max diff = %.2e",
                                   max(abs(fit_pcg$posterior_mean - fit_chol$posterior_mean))), bty = "n")

# -----------------------------------------------------------------------
# CHECK 3: Cholesky fit with explicit R_inv (diagonal weights)
# -----------------------------------------------------------------------
cat("--- Check 3: Cholesky fit with diagonal R_inv ---\n")
w       <- runif(n, 0.5, 2)
R_inv   <- Matrix::Diagonal(n, w)
fit_w   <- fit_fastblm(y, A, Q, phi = phi, R_inv = R_inv, solver = "cholesky")
cat("posterior_mean range:", range(fit_w$posterior_mean), "\n")
cat("correlation with x_true:", cor(fit_w$posterior_mean, x_true), "\n\n")

# Plot: weighted vs unweighted posterior means
plot(fit_chol$posterior_mean, fit_w$posterior_mean,
     xlab = "Unweighted posterior mean", ylab = "Weighted posterior mean",
     main = "Check 3: Effect of diagonal R_inv weighting",
     pch = 19, col = "#d01c8b", cex = 0.8)
abline(0, 1, col = "red", lty = 2)

# -----------------------------------------------------------------------
# CHECK 4: PCG fit with R_inv as operator
# -----------------------------------------------------------------------
cat("--- Check 4: PCG fit with R_inv as operator ---\n")
apply_Rinv <- function(v) w * v   # diagonal multiply
fit_pcg_w  <- fit_fastblm(y, A, Q, phi = phi,
                          R_inv = apply_Rinv, solver = "pcg")
cat("Max diff vs Cholesky weighted:", max(abs(fit_pcg_w$posterior_mean - fit_w$posterior_mean)), "\n\n")

# -----------------------------------------------------------------------
# CHECK 5: A as operator (PCG only)
# -----------------------------------------------------------------------
cat("--- Check 5: A as operator ---\n")
A_dense    <- as.matrix(A)
apply_A    <- function(v) as.numeric(A_dense %*% v)
apply_At   <- function(v) as.numeric(t(A_dense) %*% v)
fit_op     <- fit_fastblm(y, apply_A, Q, phi = phi,
                          A_t = apply_At, solver = "pcg")
cat("Max diff vs PCG matrix A:", max(abs(fit_op$posterior_mean - fit_pcg$posterior_mean)), "\n\n")

# -----------------------------------------------------------------------
# CHECK 6: posterior_se, Cholesky path
# -----------------------------------------------------------------------
cat("--- Check 6: posterior_se, Cholesky path ---\n")
se_chol <- posterior_se(fit_chol)
cat("SE range:", range(se_chol), "\n")
cat("All positive:", all(se_chol > 0), "\n\n")

# Plot: posterior mean +/- 2 SE
ord <- order(x_true)
plot(seq_len(p), fit_chol$posterior_mean[ord],
     ylim = range(c(fit_chol$posterior_mean - 2*se_chol,
                    fit_chol$posterior_mean + 2*se_chol, x_true)),
     xlab = "Coefficient index (sorted by true x)",
     ylab = "Value",
     main = "Check 6: Posterior mean ± 2 SE vs true x",
     pch = 19, col = "#2166ac", cex = 0.7)
segments(seq_len(p),
         fit_chol$posterior_mean[ord] - 2*se_chol[ord],
         seq_len(p),
         fit_chol$posterior_mean[ord] + 2*se_chol[ord],
         col = "#2166ac", lwd = 1)
points(seq_len(p), x_true[ord], pch = 4, col = "red", cex = 0.7)
legend("topleft", legend = c("Posterior mean ± 2SE", "True x"),
       col = c("#2166ac", "red"), pch = c(19, 4), bty = "n")

# -----------------------------------------------------------------------
# CHECK 7: posterior_se, Cholesky path with linear combinations
# -----------------------------------------------------------------------
cat("--- Check 7: posterior_se with A_new, Cholesky path ---\n")
A_new    <- matrix(rnorm(10 * p), 10, p)
se_lc    <- posterior_se(fit_chol, A_new = A_new)
cat("SE for 10 linear combos:", round(se_lc, 4), "\n")
cat("All positive:", all(se_lc > 0), "\n\n")

# -----------------------------------------------------------------------
# CHECK 8: posterior_se, PCG/Hutchinson path
# -----------------------------------------------------------------------
cat("--- Check 8: posterior_se, PCG/Hutchinson path ---\n")
se_pcg <- posterior_se(fit_pcg, n_probes = 200L)
cat("Correlation of PCG SE with Cholesky SE:", cor(se_pcg, se_chol), "\n")
cat("Max relative diff:", max(abs(se_pcg - se_chol) / se_chol), "\n\n")

# Plot: Hutchinson SE vs Cholesky SE
plot(se_chol, se_pcg,
     xlab = "Cholesky SE (exact)", ylab = "Hutchinson SE (stochastic)",
     main = "Check 8: Hutchinson vs Cholesky posterior SEs",
     pch = 19, col = "#e66101", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", legend = sprintf("cor = %.3f", cor(se_pcg, se_chol)), bty = "n")

# -----------------------------------------------------------------------
# CHECK 9: predict
# -----------------------------------------------------------------------
cat("--- Check 9: predict ---\n")
A_new2   <- matrix(rnorm(5 * p), 5, p)
mu_pred  <- predict(fit_chol, A_new2)
cat("Predicted means:", round(mu_pred, 4), "\n\n")

# -----------------------------------------------------------------------
# CHECK 10: REML projection
# -----------------------------------------------------------------------
cat("--- Check 10: REML projection ---\n")
X_fixed <- matrix(rnorm(n * 3), n, 3)
proj    <- reml_project(y, A, X_fixed)
cat("n_eff after projection:", proj$n_eff, "(expected", n - 3, ")\n")
cat("dim of projected A:", dim(proj$A), "\n")
cat("dim of projected y:", length(proj$y), "\n\n")

# fit on projected data
fit_reml <- fit_fastblm(proj$y, proj$A, Q, phi = phi, solver = "cholesky")
cat("REML fit posterior_mean range:", range(fit_reml$posterior_mean), "\n\n")

# -----------------------------------------------------------------------
# CHECK 11: sparse Q
# -----------------------------------------------------------------------
cat("--- Check 11: sparse banded Q (like a GMRF) ---\n")
# AR(1) precision
rho   <- 0.9
diags <- c(rep(1 + rho^2, p - 2), 1, 1) # approximate
Q_ar  <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                            diagonals = list(rep(-rho, p-1),
                                             c(1, rep(1 + rho^2, p-2), 1),
                                             rep(-rho, p-1)))
Q_ar  <- Matrix::forceSymmetric(Q_ar)
fit_ar <- fit_fastblm(y, A, Q_ar, phi = phi, solver = "cholesky")
cat("AR(1) prior fit, posterior_mean range:", range(fit_ar$posterior_mean), "\n")
cat("correlation with x_true:", cor(fit_ar$posterior_mean, x_true), "\n\n")

# Plot: AR(1) prior vs iid prior posterior means
par(mfrow = c(1, 2))
plot(x_true, fit_chol$posterior_mean,
     xlab = "True x", ylab = "Posterior mean",
     main = "iid prior",
     pch = 19, col = "#2166ac", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", legend = sprintf("cor = %.3f", cor(fit_chol$posterior_mean, x_true)), bty = "n")

plot(x_true, fit_ar$posterior_mean,
     xlab = "True x", ylab = "Posterior mean",
     main = "AR(1) prior",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", legend = sprintf("cor = %.3f", cor(fit_ar$posterior_mean, x_true)), bty = "n")
par(mfrow = c(1, 1))

cat("=== All checks complete ===\n")

# -----------------------------------------------------------------------
# CHECK 12: Woodbury fit -- p >> n case
# -----------------------------------------------------------------------
cat("--- Check 12: Woodbury fit (p >> n) ---\n")
set.seed(99)
n_w <- 100
p_w <- 500
A_w <- matrix(rnorm(n_w * p_w), n_w, p_w)
Q_w <- Matrix::Diagonal(p_w, 1)
apply_Qinv_w <- function(v) v   # Q = I so Q^{-1} = I

x_true_w <- rnorm(p_w)
y_w      <- as.numeric(A_w %*% x_true_w + rnorm(n_w))

fit_wb   <- fit_fastblm(y_w, A_w, Q_w, phi = 5, Q_inv = apply_Qinv_w, solver = "woodbury")
fit_ch_w <- fit_fastblm(y_w, A_w, Q_w, phi = 5, solver = "cholesky")
print(fit_wb)
cat("Correlation with x_true:", cor(fit_wb$posterior_mean, x_true_w), "\n")
cat("Max diff vs Cholesky:", max(abs(fit_wb$posterior_mean - fit_ch_w$posterior_mean)), "\n\n")

# -----------------------------------------------------------------------
# CHECK 13: Woodbury posterior_se vs Cholesky posterior_se
# -----------------------------------------------------------------------
cat("--- Check 13: Woodbury posterior_se vs Cholesky ---\n")
se_wb   <- posterior_se(fit_wb)
se_ch_w <- posterior_se(fit_ch_w)
cat("Correlation of SEs:", cor(se_wb, se_ch_w), "\n")
cat("Max relative diff:", max(abs(se_wb - se_ch_w) / se_ch_w), "\n")
cat("All positive:", all(se_wb > 0), "\n\n")

par(mfrow = c(1, 2))
plot(se_ch_w, se_wb,
     xlab = "Cholesky SE", ylab = "Woodbury SE",
     main = "Check 13: Woodbury vs Cholesky SE (p >> n)",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", sprintf("cor = %.4f", cor(se_wb, se_ch_w)), bty = "n")

plot(x_true_w, fit_wb$posterior_mean,
     xlab = "True x", ylab = "Posterior mean (Woodbury)",
     main = "Check 12: Woodbury fit (p >> n)",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
legend("topleft", sprintf("cor = %.3f", cor(fit_wb$posterior_mean, x_true_w)), bty = "n")
par(mfrow = c(1, 1))

dev.off()
cat("Plots saved to vignettes/sanity_checks.pdf\n")

# -----------------------------------------------------------------------
# CHECK 14: constrain -- sum to zero constraint
# -----------------------------------------------------------------------
cat("--- Check 14: constrain, sum-to-zero ---\n")
C_sum  <- matrix(1, nrow = 1, ncol = p)   # 1 x p: sum(x) = 0
fit_con <- constrain(fit_chol, C_sum)

cat("Sum of unconstrained posterior mean:", sum(fit_chol$posterior_mean), "\n")
cat("Sum of constrained posterior mean:  ", sum(fit_con$posterior_mean), "\n")
cat("(should be ~0)\n\n")

# -----------------------------------------------------------------------
# CHECK 15: constrained posterior_se should be <= unconstrained
# -----------------------------------------------------------------------
cat("--- Check 15: constrained SE <= unconstrained SE ---\n")
se_con <- posterior_se(fit_con)
cat("All constrained SE <= unconstrained SE:", all(se_con <= se_chol + 1e-10), "\n")
cat("Mean SE reduction:", mean(se_chol - se_con), "\n\n")

par(mfrow = c(1, 2))
plot(se_chol, se_con,
     xlab = "Unconstrained SE", ylab = "Constrained SE",
     main = "Check 15: Constrained vs unconstrained SE",
     pch = 19, col = "#e66101", cex = 0.8)
abline(0, 1, col = "red", lty = 2)

# -----------------------------------------------------------------------
# CHECK 16: multiple constraints
# -----------------------------------------------------------------------
cat("--- Check 16: multiple constraints ---\n")
# constrain first and last coefficient to be equal: x[1] - x[p] = 0
C_multi <- rbind(
  matrix(1, nrow = 1, ncol = p),           # sum to zero
  c(1, rep(0, p-2), -1)                     # x[1] = x[p]
)
fit_con2 <- constrain(fit_chol, C_multi)
cat("Constraint 1 (sum):", sum(fit_con2$posterior_mean), "\n")
cat("Constraint 2 (x[1]-x[p]):", fit_con2$posterior_mean[1] - fit_con2$posterior_mean[p], "\n")
cat("(both should be ~0)\n\n")

se_con2 <- posterior_se(fit_con2)
plot(se_chol, se_con2,
     xlab = "Unconstrained SE", ylab = "2-constraint SE",
     main = "Check 16: Two constraints",
     pch = 19, col = "#762a83", cex = 0.8)
abline(0, 1, col = "red", lty = 2)
par(mfrow = c(1, 1))
