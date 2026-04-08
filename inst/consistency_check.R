# Consistency check: does rho estimate converge to truth as n grows?
# Using smooth spatial A which showed identifiability
#
# This version fixes the latent simulation so that:
#   x ~ N(0, phi * sigma2e * Q^{-1})
# rather than accidentally simulating from something closer to Q^{-2}.
#
# Run with:
#   pkgload::load_all(".")
#   source("vignettes/consistency_check.R")

pkgload::load_all(".")
library(Matrix)

# -------------------------------------------------------------------
# Helper: construct 1D SAR precision
# -------------------------------------------------------------------

#' Construct a 1D SAR-induced precision matrix
#'
#' @param p Number of latent locations.
#' @param rho SAR dependence parameter.
#'
#' @return A sparse symmetric precision matrix
#'   Q = (I - rho W)'(I - rho W),
#' where W is a row-standardized nearest-neighbor matrix on a 1D chain.
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
# Helper: construct smooth observation operator
# -------------------------------------------------------------------

#' Construct a smooth spatial averaging matrix
#'
#' @param n Number of observations.
#' @param p Number of latent locations.
#' @param bw Gaussian bandwidth.
#'
#' @return An n x p row-normalized smoothing matrix.
# make_smooth_A <- function(n, p, bw = 3 / p) {
#   obs_locs <- seq(0, 1, length.out = n)
#   bas_locs <- seq(0, 1, length.out = p)
#
#   A <- matrix(0, n, p)
#
#   for (i in seq_len(n)) {
#     w <- exp(-0.5 * ((obs_locs[i] - bas_locs) / bw)^2)
#     A[i, ] <- w / sum(w)
#   }
#
#   A
# }

make_smooth_A <- function(n, p, bw = 3 / p) {


  A <- matrix(0, n, p)
  A <- A + rnorm(n*p)

  A

}

# -------------------------------------------------------------------
# Helper: draw x ~ N(0, phi * sigma2e * Q^{-1})
# -------------------------------------------------------------------

#' Draw a latent Gaussian field with covariance proportional to Q^{-1}
#'
#' @param Q Precision matrix.
#' @param phi Spatial variance multiplier.
#' @param sigma2e Observation-noise variance scale.
#'
#' @return Numeric vector x satisfying
#'   x ~ N(0, phi * sigma2e * Q^{-1}).
draw_latent_Qinv <- function(Q, phi, sigma2e) {
  R <- chol(as.matrix(Q))
  z <- rnorm(nrow(Q))
  as.numeric(backsolve(R, z)) * sqrt(phi * sigma2e)
}

# -------------------------------------------------------------------
# Model setup
# -------------------------------------------------------------------

p <- 50
rho_true <- 0.7
phi_true <- 5
sigma2e_true <- 0.5
rho_grid <- seq(0.05, 0.90, by = 0.05)
n_sizes <- c(100, 300, 1000, 3000, 10000)
n_reps <- 10L

set.seed(420)
Q_true <- make_SAR_Q(p, rho_true)
probes <- matrix(sample(c(-1L, 1L), p * 50, replace = TRUE), p, 50)

cat("=== Consistency check: smooth spatial A ===\n\n")
cat(sprintf(
  "True rho = %.1f, phi = %.1f, sigma2e = %.1f\n\n",
  rho_true, phi_true, sigma2e_true
))

pdf("vignettes/consistency_check.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# -------------------------------------------------------------------
# TEST 1: ll surface for each n -- does peak sharpen near truth?
# -------------------------------------------------------------------

cat("--- Test 1: ll surface vs n ---\n")

cols <- c("#2166ac", "#4dac26", "#e66101", "#762a83", "#d01c8b")
ll_all <- matrix(NA_real_, length(n_sizes), length(rho_grid))

for (ni in seq_along(n_sizes)) {
  n <- n_sizes[ni]
  set.seed(420)

  x_true <- draw_latent_Qinv(Q_true, phi_true, sigma2e_true)
  A <- make_smooth_A(n, p)
  y <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

  AtRinvA <- Matrix::crossprod(A)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
  AtRinvy <- as.numeric(Matrix::crossprod(A, y))
  yRinvy <- as.numeric(crossprod(y, y))

  ll_vals <- sapply(rho_grid, function(rho) {
    Q <- make_SAR_Q(p, rho)
    apply_Q <- as_apply(Q)
    logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
    prior <- list(Q = Q, Q_matrix = Q, AtRinvA_matrix = AtRinvA)

    phi_i <- fastblm:::.profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, log(0.01), log(1000),
      "cholesky", 1e-6, 4L * p, rep(0, p)
    )

    res <- fastblm:::.eval_reml_ll(
      phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L * p, rep(0, p)
    )

    res$ll
  })

  ll_all[ni, ] <- ll_vals - max(ll_vals)
  peak_rho <- rho_grid[which.max(ll_vals)]

  cat(sprintf("  n = %5d: peak rho = %.2f\n", n, peak_rho))
}

plot(
  rho_grid, ll_all[1, ],
  type = "b", pch = 19, col = cols[1],
  xlab = "rho", ylab = "REML ll (normalized)",
  main = "ll surface vs n\n(normalized to peak = 0)",
  ylim = range(ll_all, na.rm = TRUE)
)

for (ni in 2:length(n_sizes)) {
  lines(rho_grid, ll_all[ni, ], type = "b", pch = 19, col = cols[ni])
}

abline(v = rho_true, col = "red", lty = 2)
legend(
  "bottomleft",
  sprintf("n=%d", n_sizes),
  col = cols,
  pch = 19,
  lty = 1,
  bty = "n",
  cex = 0.8
)

# -------------------------------------------------------------------
# TEST 2: replicated estimation -- does variance shrink with n?
# -------------------------------------------------------------------

cat("\n--- Test 2: replicated estimation ---\n")

results <- data.frame()

for (n in n_sizes) {
  rho_ests <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100+1)

    x_rep <- draw_latent_Qinv(Q_true, phi_true, sigma2e_true)
    A <- make_smooth_A(n, p)
    y <- as.numeric(A %*% x_rep + rnorm(n, sd = sqrt(sigma2e_true)))

    AtRinvA <- Matrix::crossprod(A)
    apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
    AtRinvy <- as.numeric(Matrix::crossprod(A, y))
    yRinvy <- as.numeric(crossprod(y, y))

    ll_vals <- sapply(rho_grid, function(rho) {
      Q <- make_SAR_Q(p, rho)
      apply_Q <- as_apply(Q)
      logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
      prior <- list(Q = Q, Q_matrix = Q, AtRinvA_matrix = AtRinvA)

      phi_i <- fastblm:::.profile_phi(
        AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n, logdetQ, probes, 50L, log(0.01), log(1000),
        "cholesky", 1e-6, 4L * p, rep(0, p)
      )

      res <- fastblm:::.eval_reml_ll(
        phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L * p, rep(0, p)
      )

      res$ll
    })

    rho_ests[rep] <- rho_grid[which.max(ll_vals)]
  }

  cat(sprintf(
    "  n = %5d: mean = %.3f  sd = %.3f  bias = %.3f\n",
    n, mean(rho_ests), sd(rho_ests), mean(rho_ests) - rho_true
  ))

  results <- rbind(
    results,
    data.frame(
      n = n,
      mean = mean(rho_ests),
      sd = sd(rho_ests),
      bias = mean(rho_ests) - rho_true,
      rmse = sqrt(mean((rho_ests - rho_true)^2))
    )
  )
}

plot(
  results$n, results$mean,
  type = "b", pch = 19, col = "#2166ac", log = "x",
  xlab = "n (log scale)", ylab = "mean estimated rho",
  main = "Mean rho estimate vs n\n(should converge to true)",
  ylim = c(0.4, 1.0)
)

abline(h = rho_true, col = "red", lty = 2)
segments(
  results$n, results$mean - results$sd,
  results$n, results$mean + results$sd,
  col = "#2166ac", lwd = 2
)

legend(
  "topright",
  sprintf("True rho = %.1f", rho_true),
  col = "red",
  lty = 2,
  bty = "n"
)

plot(
  results$n, results$sd,
  type = "b", pch = 19, col = "#e66101", log = "xy",
  xlab = "n (log scale)", ylab = "SD of rho estimate (log scale)",
  main = "Estimation variance vs n\n(should shrink as 1/sqrt(n))"
)

n_ref <- results$n
lines(
  n_ref,
  results$sd[1] * sqrt(results$n[1] / n_ref),
  col = "gray",
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  c("Empirical SD", "1/sqrt(n) reference"),
  col = c("#e66101", "gray"),
  lty = c(1, 2),
  pch = c(19, NA),
  bty = "n"
)

plot(
  results$n, results$rmse,
  type = "b", pch = 19, col = "#4dac26", log = "xy",
  xlab = "n (log scale)", ylab = "RMSE (log scale)",
  main = "RMSE vs n"
)

plot(
  results$n, results$bias,
  type = "b", pch = 19, col = "#762a83", log = "x",
  xlab = "n (log scale)", ylab = "Bias (estimated - true)",
  main = "Bias vs n\n(should shrink toward 0)"
)
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))
dev.off()

cat("\nPlots saved to vignettes/consistency_check.pdf\n")
