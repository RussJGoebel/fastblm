# Spectral identifiability diagnostic for observation operators A
#
# Goal:
# Compare how different observation operators preserve latent modes.
# We use eigenvectors of a SAR precision matrix Q as a frequency-ordered basis.
#
# For each mode v_j, compute:
#   g_j = ||A v_j||^2
#
# Compare:
#   1) random A
#   2) smooth (local averaging) A
#
# Output:
#   - plots of modal gains
#   - simple summaries of identifiable cutoff

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
  as.matrix(Matrix::crossprod(S))
}

# -------------------------------------------------------------------
# A constructors
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
# Modal gains
# -------------------------------------------------------------------

compute_mode_gains <- function(A, V) {
  apply(V, 2, function(v) {
    sum((A %*% v)^2)
  })
}

# -------------------------------------------------------------------
# Summary metrics
# -------------------------------------------------------------------

summarize_gains <- function(
    gains,
    rel_threshold = 0.1,
    cumulative_threshold = 0.95
) {
  gains_rel <- gains / max(gains)
  gains_norm <- gains / sum(gains)

  last_rel_mode <- max(which(gains_rel >= rel_threshold))
  cum_cutoff <- which(cumsum(gains_norm) >= cumulative_threshold)[1]

  list(
    last_mode_above_relative_threshold = last_rel_mode,
    cumulative_energy_cutoff_mode = cum_cutoff
  )
}

# -------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------

plot_A_rows <- function(A, n_show = 6, main = "Rows of A") {
  matplot(
    t(A[seq_len(min(n_show, nrow(A))), , drop = FALSE]),
    type = "l",
    lty = 1,
    xlab = "latent index",
    ylab = "weight",
    main = main
  )
}

plot_mode_gains <- function(gains, main = "Modal gains") {
  plot(
    seq_along(gains), gains,
    type = "b",
    pch = 19,
    cex = 0.6,
    xlab = "mode index (low to high frequency)",
    ylab = expression(g[j] == group("||", A %*% v[j], "||")^2),
    main = main
  )
}

# -------------------------------------------------------------------
# Main experiment
# -------------------------------------------------------------------

p <- 100
n <- 100
rho_ref <- 0.7

set.seed(123)

Q <- make_SAR_Q(p, rho_ref)

# Eigen decomposition (frequency ordering)
eig <- eigen(Q, symmetric = TRUE)
V <- eig$vectors

# Build A matrices
A_random <- make_random_A(n, p)
A_smooth <- make_smooth_A(n, p)

# Compute gains
g_random <- compute_mode_gains(A_random, V)
g_smooth <- compute_mode_gains(A_smooth, V)

# Normalize
g_random_rel <- g_random / max(g_random)
g_smooth_rel <- g_smooth / max(g_smooth)

# Summaries
cat("=== Spectral Identifiability Summary ===\n\n")

cat("Random A:\n")
print(summarize_gains(g_random))

cat("\nSmooth A:\n")
print(summarize_gains(g_smooth))

# -------------------------------------------------------------------
# Plot results
# -------------------------------------------------------------------

pdf("test_A_spectral_modes.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# A structure
plot_A_rows(A_random, main = "Random A (rows)")
plot_A_rows(A_smooth, main = "Smooth A (rows)")

# Gains
plot_mode_gains(g_random, main = "Random A gains")
plot_mode_gains(g_smooth, main = "Smooth A gains")

# Relative comparison
plot(
  seq_len(p), g_random_rel,
  type = "b", pch = 19, cex = 0.5,
  xlab = "mode index",
  ylab = "relative gain",
  main = "Relative gains",
  ylim = range(c(g_random_rel, g_smooth_rel))
)
lines(seq_len(p), g_smooth_rel, type = "b", pch = 19, cex = 0.5, lty = 2)

legend(
  "topright",
  legend = c("Random A", "Smooth A"),
  lty = c(1, 2),
  pch = 19,
  bty = "n"
)

# Cumulative energy
plot(
  seq_len(p), cumsum(g_random / sum(g_random)),
  type = "l", lwd = 2,
  xlab = "mode index",
  ylab = "cumulative gain",
  main = "Cumulative observable energy"
)
lines(seq_len(p), cumsum(g_smooth / sum(g_smooth)), lwd = 2, lty = 2)

abline(h = 0.95, col = "gray", lty = 3)

legend(
  "bottomright",
  legend = c("Random A", "Smooth A", "95% cutoff"),
  lty = c(1, 2, 3),
  lwd = c(2, 2, 1),
  bty = "n"
)

par(mfrow = c(1, 1))
dev.off()

cat("\nSaved plots to test_A_spectral_modes.pdf\n")
