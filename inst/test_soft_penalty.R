# =============================================================================
# test_soft_penalty.R
#
# Tests that the hard RSR constraint is the limit of a soft quadratic penalty.
#
# The hard constraint enforced by constrain(fit, C) is:
#   C * x = 0  (via Schur complement correction)
#
# This is equivalent to the limit as lambda -> Inf of adding a soft penalty:
#   Q_penalised = Q + lambda * t(C) %*% C
#
# As lambda increases, fit_fastblm(Q_penalised) should converge to
# constrain(fit, C) in the posterior mean.
#
# Usage:
#   source("utils.R"); source("solvers.R"); source("fit.R"); source("constrain.R")
#   source("test_soft_penalty.R")
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
# Setup
# =============================================================================
cat("=====================================================\n")
cat(" Soft penalty -> hard constraint convergence test\n")
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

# covariates
set.seed(1)
cov_latent <- matrix(rnorm(p), m, m)
for (i in 1:3)
  cov_latent <- (cov_latent +
                   rbind(cov_latent[-1,], cov_latent[nrow(cov_latent),]) +
                   rbind(cov_latent[1,],  cov_latent[-nrow(cov_latent),])) / 3
X_latent <- cbind(1, scale(as.vector(cov_latent)))
X_obs    <- as.matrix(A %*% X_latent)
q        <- ncol(X_obs)

# simulate data
r <- simulate_sar(W, rho, phi, sigma2e)
x <- as.numeric(X_latent %*% beta) + r
y <- as.numeric(A %*% x) + rnorm(n, sd = sqrt(sigma2e))

# block design
A_full  <- cbind(X_obs, as.matrix(A))
epsilon <- 1e-6
S       <- Matrix::Diagonal(p) - rho * W
Q_sar   <- Matrix::forceSymmetric(Matrix::crossprod(S))
Q_full  <- Matrix::bdiag(Matrix::Diagonal(q, epsilon), Q_sar)

# constraint: C_full = [0 | t(X_obs) %*% A]
C_r    <- as.matrix(t(X_obs) %*% A)
C_full <- cbind(matrix(0, q, q), C_r)

cat(sprintf("Grid: %dx%d = %d pixels,  n_obs = %d\n", m, m, p, n))
cat(sprintf("Block design: %d x %d,  constraint: %d x %d\n\n",
            nrow(A_full), ncol(A_full), nrow(C_full), ncol(C_full)))

# =============================================================================
# Hard constraint reference
# =============================================================================
fit_unc  <- fit_fastblm(y, A_full, Q_full, phi = phi, solver = "cholesky")
fit_hard <- constrain(fit_unc, C_full)
mu_hard  <- fit_hard$posterior_mean

cat(sprintf("Hard constraint: max|C * mu| = %.2e  (should be ~0)\n\n",
            max(abs(C_full %*% mu_hard))))

# =============================================================================
# Soft penalty sweep
# =============================================================================
CtC     <- Matrix::crossprod(C_full)
lambdas <- c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e8)
diffs   <- numeric(length(lambdas))
viols   <- numeric(length(lambdas))

cat(sprintf("  %-10s  %-24s  %-20s\n",
            "lambda", "max|mu_soft - mu_hard|", "max|C * mu_soft|"))
cat(sprintf("  %s\n", strrep("-", 58)))

for (i in seq_along(lambdas)) {
  lam      <- lambdas[i]
  Q_soft   <- Q_full + lam * CtC
  fit_soft <- fit_fastblm(y, A_full, Q_soft, phi = phi, solver = "cholesky")
  diffs[i] <- max(abs(fit_soft$posterior_mean - mu_hard))
  viols[i] <- max(abs(C_full %*% fit_soft$posterior_mean))
  cat(sprintf("  %-10.0e  %-24.8f  %-20.8f\n", lam, diffs[i], viols[i]))
}

# =============================================================================
# Pass criteria
# =============================================================================
cat("\n")
monotone <- all(diff(diffs) < 0)
converged <- diffs[length(diffs)] < 1e-3

cat(sprintf("  Differences monotonically decreasing: %s\n", monotone))
cat(sprintf("  Final error < 1e-3: %s  (actual: %.2e)\n", converged, diffs[length(diffs)]))

pass <- monotone && converged
cat(sprintf("\nOverall: %s\n", if(pass) "PASS" else "FAIL"))
