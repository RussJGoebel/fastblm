# Profile Fisher information for rho under different observation operators A
#
# Goal:
# For the covariance model
#   Sigma(rho, phi) = I + phi * A Q(rho)^{-1} A'
# compute the Fisher blocks
#   I_rhorho, I_rhophi, I_phiphi
# and then the profile Fisher information
#   I_profile(rho) = I_rhorho - I_rhophi^2 / I_phiphi
#
# Compare:
#   1) random A
#   2) smooth A
#
# This is a simulation-free identifiability diagnostic for rho after
# profiling out phi.
#
# Notes:
# - sigma2_e is omitted because Fisher info for rho is invariant to an
#   overall covariance scale multiplier.
# - Sigma_rho is approximated by finite differences on the rho grid.

library(Matrix)

# -------------------------------------------------------------------
# SAR precision
# -------------------------------------------------------------------

make_SAR_Q <- function(p, rho) {
  i_idx <- c(2:p, 1:(p - 1))
  j_idx <- c(1:(p - 1), 2:p)
  row_sums <- c(1, rep(2, p - 2), 1)
  x_vals <- 1 / row_sums[i_idx]

  W <- Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_vals,
    dims = c(p, p)
  )

  S <- Matrix::Diagonal(p) - rho * W
  Matrix::crossprod(S)
}

# -------------------------------------------------------------------
# Observation operators
# -------------------------------------------------------------------

make_random_A <- function(n, p) {
  A <- matrix(rnorm(n * p), n, p)
  row_norms <- sqrt(rowSums(A^2))
  A / row_norms
}

make_smooth_A <- function(n, p, bw = 3 / p) {
  obs_locs <- seq(0, 1, length.out = n)
  bas_locs <- seq(0, 1, length.out = p)

  A <- matrix(0, n, p)

  for (i in seq_len(n)) {
    w <- exp(-0.5 * ((obs_locs[i] - bas_locs) / bw)^2)
    A[i, ] <- w / sum(w)
  }

  A
}

# -------------------------------------------------------------------
# Compute AQ^{-1}A' without forming Q^{-1}
# -------------------------------------------------------------------

compute_sigma_obs <- function(A, Q) {
  X <- solve(as.matrix(Q), t(A))
  A %*% X
}

# -------------------------------------------------------------------
# Finite-difference derivative of a list of matrices
# -------------------------------------------------------------------

compute_matrix_derivative <- function(M_list, x_grid) {
  m <- length(M_list)
  dM <- vector("list", m)

  for (i in 2:(m - 1)) {
    dx <- x_grid[i + 1] - x_grid[i - 1]
    dM[[i]] <- (M_list[[i + 1]] - M_list[[i - 1]]) / dx
  }

  dM[[1]] <- (M_list[[2]] - M_list[[1]]) / (x_grid[2] - x_grid[1])
  dM[[m]] <- (M_list[[m]] - M_list[[m - 1]]) / (x_grid[m] - x_grid[m - 1])

  dM
}

# -------------------------------------------------------------------
# Fisher blocks for one rho
# -------------------------------------------------------------------

compute_fisher_blocks <- function(Sigma, Sigma_rho, Sigma_phi) {
  Sigma_inv <- solve(Sigma)

  M_rho <- Sigma_inv %*% Sigma_rho
  M_phi <- Sigma_inv %*% Sigma_phi

  I_rhorho <- 0.5 * sum(M_rho * t(M_rho))
  I_rhophi <- 0.5 * sum(M_rho * t(M_phi))
  I_phiphi <- 0.5 * sum(M_phi * t(M_phi))

  c(
    I_rhorho = I_rhorho,
    I_rhophi = I_rhophi,
    I_phiphi = I_phiphi
  )
}

# -------------------------------------------------------------------
# Compute profile Fisher curve for a given A
# -------------------------------------------------------------------

compute_profile_fisher_curve <- function(A, rho_grid, phi) {
  # Step 1: compute Sigma_obs(rho) = A Q(rho)^{-1} A'
  Sigma_obs_list <- lapply(rho_grid, function(rho) {
    Q <- make_SAR_Q(ncol(A), rho)
    compute_sigma_obs(A, Q)
  })

  # Step 2: derivative wrt rho of Sigma_obs(rho)
  dSigma_obs_list <- compute_matrix_derivative(Sigma_obs_list, rho_grid)

  # Step 3: full covariance and derivatives
  # Sigma(rho, phi) = I + phi * Sigma_obs(rho)
  m <- length(rho_grid)

  I_rhorho <- numeric(m)
  I_rhophi <- numeric(m)
  I_phiphi <- numeric(m)
  I_profile <- numeric(m)

  for (i in seq_len(m)) {
    Sigma_obs <- Sigma_obs_list[[i]]
    dSigma_obs <- dSigma_obs_list[[i]]

    Sigma <- diag(nrow(A)) + phi * Sigma_obs
    Sigma_rho <- phi * dSigma_obs
    Sigma_phi <- Sigma_obs

    blocks <- compute_fisher_blocks(Sigma, Sigma_rho, Sigma_phi)

    I_rhorho[i] <- blocks["I_rhorho"]
    I_rhophi[i] <- blocks["I_rhophi"]
    I_phiphi[i] <- blocks["I_phiphi"]

    I_profile[i] <- I_rhorho[i] - I_rhophi[i]^2 / I_phiphi[i]
  }

  list(
    rho_grid = rho_grid,
    I_rhorho = I_rhorho,
    I_rhophi = I_rhophi,
    I_phiphi = I_phiphi,
    I_profile = I_profile,
    Sigma_obs_list = Sigma_obs_list
  )
}

# -------------------------------------------------------------------
# Main setup
# -------------------------------------------------------------------

set.seed(123)

p <- 100
n <- 100
phi <- 5
rho_grid <- seq(0.05, 0.95, by = 0.05)

A_random <- make_random_A(n, p)
A_smooth <- make_smooth_A(n, p, bw = 3 / p)

# -------------------------------------------------------------------
# Compute profile Fisher curves
# -------------------------------------------------------------------

cat("Computing profile Fisher for random A...\n")
fish_random <- compute_profile_fisher_curve(A_random, rho_grid, phi)

cat("Computing profile Fisher for smooth A...\n")
fish_smooth <- compute_profile_fisher_curve(A_smooth, rho_grid, phi)

# Relative versions for shape comparison
rel <- function(x) {
  x / max(x, na.rm = TRUE)
}

Iprof_random_rel <- rel(fish_random$I_profile)
Iprof_smooth_rel <- rel(fish_smooth$I_profile)

Irhorho_random_rel <- rel(fish_random$I_rhorho)
Irhorho_smooth_rel <- rel(fish_smooth$I_rhorho)

# -------------------------------------------------------------------
# Print summaries
# -------------------------------------------------------------------

cat("\n=== Profile Fisher summaries ===\n")
cat(sprintf("Random A: max I_profile = %.6f\n", max(fish_random$I_profile)))
cat(sprintf("Smooth A: max I_profile = %.6f\n", max(fish_smooth$I_profile)))

cat(sprintf(
  "Random A: min/max I_profile = [%.6f, %.6f]\n",
  min(fish_random$I_profile), max(fish_random$I_profile)
))
cat(sprintf(
  "Smooth A: min/max I_profile = [%.6f, %.6f]\n",
  min(fish_smooth$I_profile), max(fish_smooth$I_profile)
))

# -------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------

pdf("profile_fisher_info.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# Show A rows
matplot(
  t(A_random[1:6, , drop = FALSE]),
  type = "l",
  lty = 1,
  xlab = "latent index",
  ylab = "weight",
  main = "Random A (rows)"
)

matplot(
  t(A_smooth[1:6, , drop = FALSE]),
  type = "l",
  lty = 1,
  xlab = "latent index",
  ylab = "weight",
  main = "Smooth A (rows)"
)

# Raw I_rhorho
plot(
  rho_grid, fish_random$I_rhorho,
  type = "b", pch = 19,
  xlab = expression(rho),
  ylab = expression(I[rho * rho]),
  main = expression("Raw Fisher block " * I[rho * rho]),
  ylim = range(c(fish_random$I_rhorho, fish_smooth$I_rhorho))
)
lines(rho_grid, fish_smooth$I_rhorho, type = "b", pch = 19, lty = 2)
legend(
  "topleft",
  legend = c("Random A", "Smooth A"),
  lty = c(1, 2),
  pch = 19,
  bty = "n"
)

# Raw profile Fisher
plot(
  rho_grid, fish_random$I_profile,
  type = "b", pch = 19,
  xlab = expression(rho),
  ylab = expression(I[profile](rho)),
  main = "Profile Fisher information",
  ylim = range(c(fish_random$I_profile, fish_smooth$I_profile))
)
lines(rho_grid, fish_smooth$I_profile, type = "b", pch = 19, lty = 2)
legend(
  "topleft",
  legend = c("Random A", "Smooth A"),
  lty = c(1, 2),
  pch = 19,
  bty = "n"
)

# Relative I_rhorho
plot(
  rho_grid, Irhorho_random_rel,
  type = "b", pch = 19,
  xlab = expression(rho),
  ylab = "relative value",
  main = expression("Relative " * I[rho * rho]),
  ylim = range(c(Irhorho_random_rel, Irhorho_smooth_rel))
)
lines(rho_grid, Irhorho_smooth_rel, type = "b", pch = 19, lty = 2)
legend(
  "topleft",
  legend = c("Random A", "Smooth A"),
  lty = c(1, 2),
  pch = 19,
  bty = "n"
)

# Relative profile Fisher
plot(
  rho_grid, Iprof_random_rel,
  type = "b", pch = 19,
  xlab = expression(rho),
  ylab = "relative value",
  main = "Relative profile Fisher",
  ylim = range(c(Iprof_random_rel, Iprof_smooth_rel))
)
lines(rho_grid, Iprof_smooth_rel, type = "b", pch = 19, lty = 2)
legend(
  "topleft",
  legend = c("Random A", "Smooth A"),
  lty = c(1, 2),
  pch = 19,
  bty = "n"
)

par(mfrow = c(1, 1))
dev.off()

cat("\nSaved plots to profile_fisher_info.pdf\n")
