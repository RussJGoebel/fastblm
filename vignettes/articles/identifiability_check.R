# Identifiability check for SAR rho
# Key question: does more data help identify rho?
# If the ll surface stays monotonic regardless of n, it's a model issue.
# If it becomes unimodal with large n, it's a finite sample issue.
# Run with: pkgload::load_all("."); source("vignettes/identifiability_check.R")

pkgload::load_all(".")
library(Matrix)

make_SAR_Q <- function(p, rho) {
  # W: row-normalized adjacency for 1D chain
  # boundary nodes have 1 neighbor, interior nodes have 2
  i_idx    <- c(2:p, 1:(p-1))
  j_idx    <- c(1:(p-1), 2:p)
  row_sums <- c(1, rep(2, p-2), 1)
  x_vals   <- 1 / row_sums[i_idx]
  W <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_vals, dims = c(p, p))

  # Q = (I - rho W)' (I - rho W)
  S <- Matrix::Diagonal(p) - rho * W
  Matrix::forceSymmetric(Matrix::crossprod(S))
}

p         <- 100
rho_true  <- 0.7
phi_true  <- 8
sigma2e_true <- 0.5
rho_grid  <- seq(0.05, 0.95, by = 0.05)

set.seed(42)
Q_true <- make_SAR_Q(p, rho_true)
x_true <- as.numeric(Matrix::solve(Matrix::Cholesky(Q_true), rnorm(p))) *
  sqrt(phi_true * sigma2e_true)

pdf("vignettes/identifiability_check.pdf", width = 14, height = 14)
par(mfrow = c(2, 3))

n_sizes <- c(100, 300, 1000, 3000, 10000)

for (n in n_sizes) {
  cat(sprintf("n = %d\n", n))
  set.seed(42)
  A <- matrix(rnorm(n * p), n, p) / sqrt(p)
  y <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

  Rinv          <- resolve_Rinv(NULL, n)
  AtRinvy       <- as.numeric(Matrix::crossprod(A, y))
  yRinvy        <- as.numeric(crossprod(y, y))
  probes        <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)
  AtRinvA       <- Matrix::crossprod(A)   # precompute once
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)

  ref_ll <- sapply(rho_grid, function(rho) {
    Q       <- make_SAR_Q(p, rho)
    apply_Q <- as_apply(Q)
    logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
    prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA)

    phi_i <- fastblm:::.profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, log(0.01), log(1000),
      "cholesky", 1e-6, 4L*p, rep(0,p)
    )

    res <- fastblm:::.eval_reml_ll(
      phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res$ll
  })

  plot(rho_grid, ref_ll,
       type = "b", pch = 19, col = "#2166ac",
       xlab = "rho", ylab = "REML ll",
       main = sprintf("n = %d", n))
  abline(v = rho_true, col = "red", lty = 2)
  legend("bottomright", sprintf("True rho = %.1f", rho_true),
         col = "red", lty = 2, bty = "n")
}

# -----------------------------------------------------------------------
# Also try: spatial design where A has local structure
# Each observation is a local average of nearby basis functions
# This should make rho more identifiable
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Spatial design with varying n
# -----------------------------------------------------------------------
cat("\nTrying spatial design matrix with varying n...\n")

n_sizes_spatial <- c(300, 1000, 3000, 10000)

par(mfrow = c(2, 2))

for (n in n_sizes_spatial) {
  cat(sprintf("  n = %d\n", n))
  set.seed(42)
  rows <- rep(seq_len(n), each = 5)
  cols <- unlist(lapply(seq_len(n), function(i) {
    start <- sample(seq_len(p - 4), 1)
    start:(start + 4)
  }))
  vals <- rep(1/5, n * 5)
  A_spatial <- as.matrix(Matrix::sparseMatrix(i=rows, j=cols, x=vals, dims=c(n, p)))
  y_spatial  <- as.numeric(A_spatial %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

  Rinv          <- resolve_Rinv(NULL, n)
  AtRinvy       <- as.numeric(Matrix::crossprod(A_spatial, y_spatial))
  yRinvy        <- as.numeric(crossprod(y_spatial, y_spatial))
  probes        <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)
  AtRinvA       <- Matrix::crossprod(A_spatial)   # precompute once
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)

  ll_spatial <- sapply(rho_grid, function(rho) {
    Q       <- make_SAR_Q(p, rho)
    apply_Q <- as_apply(Q)
    logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
    prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA)

    phi_i <- fastblm:::.profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, log(0.01), log(1000),
      "cholesky", 1e-6, 4L*p, rep(0,p)
    )

    res <- fastblm:::.eval_reml_ll(
      phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res$ll
  })

  peak_rho <- rho_grid[which.max(ll_spatial)]
  plot(rho_grid, ll_spatial,
       type = "b", pch = 19, col = "#e66101",
       xlab = "rho", ylab = "REML ll",
       main = sprintf("Spatial design, n = %d\npeak at rho = %.2f", n, peak_rho))
  abline(v = rho_true,  col = "red",      lty = 2)
  abline(v = peak_rho,  col = "#4dac26",  lty = 2)
  legend("bottomright",
         legend = c(sprintf("True rho = %.1f", rho_true),
                    sprintf("Peak rho = %.2f", peak_rho)),
         col = c("red", "#4dac26"), lty = 2, bty = "n")
}

par(mfrow = c(1,1))

par(mfrow = c(1,1))
dev.off()
cat("Plots saved to vignettes/identifiability_check.pdf\n")
