pkgload::load_all(".")
# =============================================================================
# test_constraint.R
#
# Tests constraint_matrix / X_cov arguments in tune_reml() and tune_cv().
#
# Full model:
#   y = X_obs*beta + A*r + eps,   eps ~ N(0, sigma2e * I)
#   r ~ SAR(rho, phi)
#
# Design matrix: A_full = [X_obs | A_spatial]
# Precision:     Q_full = blockdiag(epsilon*I_q, Q_sar)  (flat prior on beta)
#
# RSR constraint on [beta; r]:
#   C_full = [0_{q x q} | t(X_obs) %*% A]
# Forces A*r orthogonal to X_obs, making beta identifiable.
#
# Tests:
#   1. MECHANICAL:   C_full * posterior_mean ~ 0 after constraining
#   2. PHI RECOVERY: constrained REML recovers phi near truth
#   3. CONSTRAINT EFFECT: constrained and unconstrained phi differ
#   4. REML vs CV:   constrained REML and CV agree on phi
#   5. BETA RECOVERY: constrained model recovers covariate beta better
#
# Usage:
#   source("utils.R"); source("solvers.R"); source("fit.R")
#   source("reml.R"); source("cv.R"); source("constrain.R")
#   source("test_constraint.R")
# =============================================================================

set.seed(42)

# ---- helpers ----------------------------------------------------------------

make_grid_W <- function(m) {
  n  <- m * m
  rc <- which(matrix(TRUE, m, m), arr.ind = TRUE)
  rows <- integer(0); cols <- integer(0)
  for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next
    ni <- rc[, 1] + di;  nj <- rc[, 2] + dj
    ok <- ni >= 1 & ni <= m & nj >= 1 & nj <= m
    rows <- c(rows, which(ok))
    cols <- c(cols, (ni[ok] - 1L) * m + nj[ok])
  }
  W_raw <- Matrix::sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, n))
  Matrix::Diagonal(x = 1 / Matrix::rowSums(W_raw)) %*% W_raw
}

simulate_sar <- function(W, rho, phi, sigma2e) {
  p <- nrow(W)
  M <- Matrix::Diagonal(p) - rho * W
  as.numeric(Matrix::solve(M, rnorm(p))) * sqrt(phi * sigma2e)
}

make_agg_matrix <- function(p, n_obs, fsize) {
  stride <- max(1L, floor(p / n_obs))
  starts <- seq(1L, by = stride, length.out = n_obs)
  ends   <- pmin(starts + fsize - 1L, p)
  rows   <- rep(seq_len(n_obs), times = ends - starts + 1L)
  cols   <- unlist(mapply(seq, starts, ends, SIMPLIFY = FALSE))
  vals   <- rep(1 / (ends - starts + 1L), times = ends - starts + 1L)
  Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_obs, p))
}

# =============================================================================
# Shared setup
# =============================================================================
cat("=====================================================\n")
cat(" Constraint test: tune_reml + tune_cv\n")
cat("=====================================================\n\n")

m       <- 20
p       <- m * m
rho     <- 0.99
phi     <- 5.0
sigma2e <- 1.0
beta    <- c(2.0, -1.5)

W <- make_grid_W(m)
A <- make_agg_matrix(p, n_obs = 200L, fsize = 8L)
n <- nrow(A)

cat(sprintf("Grid %dx%d = %d pixels,  n_obs = %d\n", m, m, p, n))
cat(sprintf("True: rho=%.2f  phi=%.1f  sigma2e=%.1f\n\n", rho, phi, sigma2e))

# --- covariates --------------------------------------------------------------
set.seed(1)
cov_latent <- matrix(rnorm(p), m, m)
for (i in 1:3)
  cov_latent <- (cov_latent +
                   rbind(cov_latent[-1,], cov_latent[nrow(cov_latent),]) +
                   rbind(cov_latent[1,],  cov_latent[-nrow(cov_latent),])) / 3
X_latent <- cbind(1, scale(as.vector(cov_latent)))   # p x 2
X_obs    <- as.matrix(A %*% X_latent)                # n x 2
q        <- ncol(X_obs)

# --- simulate data -----------------------------------------------------------
r <- simulate_sar(W, rho, phi, sigma2e)
x <- as.numeric(X_latent %*% beta) + r
y <- as.numeric(A %*% x) + rnorm(n, sd = sqrt(sigma2e))

# --- block design and precision ----------------------------------------------
A_full  <- cbind(X_obs, as.matrix(A))
epsilon <- 1e-6
S       <- Matrix::Diagonal(p) - rho * W
Q_sar   <- Matrix::forceSymmetric(Matrix::crossprod(S))
Q_full  <- Matrix::bdiag(Matrix::Diagonal(q, epsilon), Q_sar)
ld_Q    <- 2 * as.numeric(Matrix::determinant(S, logarithm = TRUE)$modulus)

Q_fun_full <- function(theta) {
  list(Q = Q_full, Q_matrix = Q_full, log_det_Q = ld_Q)
}

# --- constraint --------------------------------------------------------------
C_r    <- as.matrix(t(X_obs) %*% A)
C_full <- cbind(matrix(0, q, q), C_r)

cat(sprintf("Block design: A_full is %d x %d\n", nrow(A_full), ncol(A_full)))
cat(sprintf("Constraint C_full: %d x %d\n\n", nrow(C_full), ncol(C_full)))

# =============================================================================
# Test 1: MECHANICAL
# =============================================================================
cat("--- Test 1: mechanical constraint satisfaction ---\n\n")

fit_unc <- fit_fastblm(y, A_full, Q_full, phi = phi, solver = "cholesky")
fit_con <- constrain(fit_unc, C_full)

resid_before <- C_full %*% fit_unc$posterior_mean
resid_after  <- C_full %*% fit_con$posterior_mean

cat(sprintf("  |C * mu| before: %.2e  %.2e\n",
            abs(resid_before[1]), abs(resid_before[2])))
cat(sprintf("  |C * mu| after : %.2e  %.2e\n",
            abs(resid_after[1]),  abs(resid_after[2])))

tol   <- 1e-6
pass1 <- all(abs(resid_after) < tol)
cat(sprintf("  PASS: %s  (tol = %.0e)\n\n", pass1, tol))

# =============================================================================
# Test 2: PHI RECOVERY
# =============================================================================
cat("--- Test 2: constrained REML phi recovery ---\n\n")

tuned_reml_con <- tune_reml(
  y                 = y,
  A                 = A_full,
  Q_fun             = Q_fun_full,
  theta_init        = numeric(0),
  logdet_method     = "cholesky",
  constraint_matrix = C_full,
  verbose           = FALSE
)

cat(sprintf("  True phi : %.4f  (tau = %.4f)\n", phi, 1/phi))
cat(sprintf("  REML phi : %.4f  (tau = %.4f)\n",
            tuned_reml_con$phi, 1/tuned_reml_con$phi))
phi_err <- abs(tuned_reml_con$phi - phi) / phi
cat(sprintf("  Relative error: %.1f%%\n", 100 * phi_err))
pass2 <- phi_err < 0.3
cat(sprintf("  PASS: %s  (< 30%% relative error)\n\n", pass2))

# =============================================================================
# Test 3: CONSTRAINT EFFECT
# =============================================================================
cat("--- Test 3: constraint changes phi estimate ---\n\n")

tuned_reml_nocon <- tune_reml(
  y             = y,
  A             = A_full,
  Q_fun         = Q_fun_full,
  theta_init    = numeric(0),
  logdet_method = "cholesky",
  verbose       = FALSE
)

cat(sprintf("  REML phi unconstrained: %.4f  (tau = %.4f)\n",
            tuned_reml_nocon$phi, 1/tuned_reml_nocon$phi))
cat(sprintf("  REML phi constrained  : %.4f  (tau = %.4f)\n",
            tuned_reml_con$phi,   1/tuned_reml_con$phi))
pass3 <- abs(tuned_reml_con$phi - tuned_reml_nocon$phi) > 0.01
cat(sprintf("  PASS (estimates differ): %s\n\n", pass3))

# =============================================================================
# Test 4: REML vs CV
# =============================================================================
cat("--- Test 4: constrained REML vs CV agree on phi ---\n\n")

tuned_cv_con <- tune_cv(
  y          = y,
  A          = A_full,
  Q_fun      = Q_fun_full,
  theta_init = numeric(0),
  k          = 10L,
  X_cov      = X_obs,
  seed       = 1L,
  verbose    = FALSE
)

cat(sprintf("  REML phi : %.4f  (tau = %.4f)\n",
            tuned_reml_con$phi, 1/tuned_reml_con$phi))
cat(sprintf("  CV   phi : %.4f  (tau = %.4f)\n",
            tuned_cv_con$phi,   1/tuned_cv_con$phi))
ratio <- tuned_reml_con$phi / tuned_cv_con$phi
cat(sprintf("  Ratio REML/CV: %.3f\n", ratio))
pass4 <- ratio > 0.5 && ratio < 2.0
cat(sprintf("  PASS (within 2x): %s\n\n", pass4))

# =============================================================================
# Test 5: BETA RECOVERY
#
# Large covariate effect (beta = c(2, 8)).
# Unconstrained: A*r absorbs covariate signal -> beta biased toward zero.
# Constrained:   A*r forced orthogonal to X_obs -> beta correctly identified.
# =============================================================================
cat("--- Test 5: beta recovery (RSR fixes spatial confounding) ---\n\n")

beta_strong <- c(2.0, 8.0)

set.seed(99)
r_s_raw <- simulate_sar(W, rho, phi, sigma2e)

# Project r_s onto the constraint-satisfying subspace:
# remove the component of A*r that lies in the column space of X_obs.
# This gives the "true" r as seen by the constrained model.
Ar_raw   <- as.numeric(A %*% r_s_raw)
H_X      <- X_obs %*% solve(t(X_obs) %*% X_obs, t(X_obs))   # hat matrix for X_obs
r_s      <- r_s_raw - as.numeric(Matrix::solve(
  Matrix::crossprod(A) + Matrix::Diagonal(p, 1e-8),
  Matrix::crossprod(A, H_X %*% Ar_raw)
))   # project out X_obs component from r via A

x_s <- as.numeric(X_latent %*% beta_strong) + r_s
y_s <- as.numeric(A %*% x_s) + rnorm(n, sd = sqrt(sigma2e))

tuned_s <- tune_reml(
  y                 = y_s,
  A                 = A_full,
  Q_fun             = Q_fun_full,
  theta_init        = numeric(0),
  logdet_method     = "cholesky",
  constraint_matrix = C_full,
  verbose           = FALSE
)

fit_nocon_s <- fit_fastblm(y_s, A_full, Q_full,
                           phi = tuned_s$phi, solver = "cholesky")
fit_con_s   <- constrain(fit_nocon_s, C_full)

beta_nocon <- fit_nocon_s$posterior_mean[seq_len(q)]
beta_con   <- fit_con_s$posterior_mean[seq_len(q)]

cat(sprintf("  True beta          : intercept=%.3f  covariate=%.3f\n",
            beta_strong[1], beta_strong[2]))
cat(sprintf("  beta unconstrained : intercept=%.3f  covariate=%.3f\n",
            beta_nocon[1], beta_nocon[2]))
cat(sprintf("  beta constrained   : intercept=%.3f  covariate=%.3f\n",
            beta_con[1], beta_con[2]))

err_cov_nocon <- abs(beta_nocon[2] - beta_strong[2])
err_cov_con   <- abs(beta_con[2]   - beta_strong[2])

cat(sprintf("\n  |error| covariate unconstrained: %.4f\n", err_cov_nocon))
cat(sprintf("  |error| covariate constrained  : %.4f\n", err_cov_con))
cat(sprintf("  Improvement: %.1f%%\n",
            100 * (err_cov_nocon - err_cov_con) / err_cov_nocon))

pass5 <- err_cov_con < err_cov_nocon
cat(sprintf("  PASS (constrained recovers covariate better): %s\n\n", pass5))

# =============================================================================
# Summary
# =============================================================================
cat("=====================================================\n")
cat(" SUMMARY\n")
cat("=====================================================\n")
cat(sprintf("  Test 1 (constraint satisfied)   : %s\n", if(pass1) "PASS" else "FAIL"))
cat(sprintf("  Test 2 (phi recovery <30%%)      : %s\n", if(pass2) "PASS" else "FAIL"))
cat(sprintf("  Test 3 (constraint changes phi) : %s\n", if(pass3) "PASS" else "FAIL"))
cat(sprintf("  Test 4 (REML/CV agree <2x)      : %s\n", if(pass4) "PASS" else "FAIL"))
cat(sprintf("  Test 5 (beta recovery improved) : %s\n", if(pass5) "PASS" else "FAIL"))
cat("\n")
all_pass <- pass1 && pass2 && pass3 && pass4 && pass5
cat(sprintf("Overall: %s\n", if(all_pass) "ALL PASS" else "SOME FAILURES"))
