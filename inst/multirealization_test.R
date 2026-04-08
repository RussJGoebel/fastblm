# Can we identify rho from multiple independent field realizations?
# Key idea:
# With K independent realizations y_1, ..., y_K from the same covariance model,
# the empirical covariance of y should identify the shape of
#   A Q(rho)^{-1} A'
# if the simulation and fitting model actually match.
#
# This version fixes the latent sampler so that:
#   x ~ N(0, phi * sigma2e * Q^{-1})
# rather than accidentally producing something closer to Q^{-2}.
#
# Run with:
#   pkgload::load_all(".")
#   source("vignettes/multi_realization_test.R")

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
# Helper: smooth observation operator
# -------------------------------------------------------------------

#' Construct a smooth observation matrix A
#'
#' @param n Number of observations per realization.
#' @param p Number of latent locations.
#' @param bw Kernel bandwidth.
#'
#' @return An n x p matrix whose rows average nearby latent positions.
make_smooth_A <- function(n, p, bw = 3 / p) {
  obs <- seq(0, 1, length.out = n)
  bas <- seq(0, 1, length.out = p)

  D <- outer(
    obs,
    bas,
    function(o, b) exp(-0.5 * ((o - b) / bw)^2)
  )

  D / rowSums(D)
}

# -------------------------------------------------------------------
# Helper: latent Gaussian sampler with covariance Q^{-1}
# -------------------------------------------------------------------

#' Draw a latent field with covariance proportional to Q^{-1}
#'
#' @param Q Precision matrix.
#' @param phi Spatial variance multiplier.
#' @param sigma2e Noise variance scale.
#'
#' @return Numeric vector x with
#'   x ~ N(0, phi * sigma2e * Q^{-1}).
draw_latent_Qinv <- function(Q, phi, sigma2e) {
  R <- chol(as.matrix(Q))
  z <- rnorm(nrow(Q))
  as.numeric(backsolve(R, z)) * sqrt(phi * sigma2e)
}

# -------------------------------------------------------------------
# Helper: simulate K realizations of y
# -------------------------------------------------------------------

#' Simulate multiple independent realizations
#'
#' @param K Number of realizations.
#' @param Q_true True latent precision.
#' @param A Observation matrix.
#' @param phi_true True phi.
#' @param sigma2e_true True noise variance.
#'
#' @return n x K matrix of observed realizations.
simulate_Y_mat <- function(K, Q_true, A, phi_true, sigma2e_true) {
  n <- nrow(A)

  Y_mat <- sapply(seq_len(K), function(k) {
    x_k <- draw_latent_Qinv(Q_true, phi_true, sigma2e_true)
    as.numeric(A %*% x_k + rnorm(n, sd = sqrt(sigma2e_true)))
  })

  if (K == 1) {
    Y_mat <- matrix(Y_mat, ncol = 1)
  }

  Y_mat
}

# -------------------------------------------------------------------
# Model settings
# -------------------------------------------------------------------

p <- 50
rho_true <- 0.7
phi_true <- 5
sigma2e_true <- 0.5
n_per_field <- 50
rho_grid <- seq(0.05, 0.90, by = 0.025)

Q_true <- make_SAR_Q(p, rho_true)
A <- make_smooth_A(n_per_field, p)

# -------------------------------------------------------------------
# Precompute Q quantities on the rho grid
# -------------------------------------------------------------------

Q_precomp <- lapply(rho_grid, function(rho) {
  Q <- make_SAR_Q(p, rho)
  Qinv <- solve(as.matrix(Q))
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)

  list(
    Q = Q,
    Qinv = Qinv,
    logdetQ = logdetQ
  )
})

# -------------------------------------------------------------------
# Marginal log-likelihood with phi fixed
# y_k ~ N(0, sigma2e * [I + phi A Q^{-1} A'])
# Sigma2e is profiled out jointly across all K realizations.
# -------------------------------------------------------------------

#' Marginal log-likelihood across K independent realizations
#'
#' @param rho_idx Index into rho_grid / Q_precomp.
#' @param Y_mat n x K matrix of observations.
#' @param A Observation matrix.
#' @param phi Spatial variance multiplier.
#' @param n Number of observations per realization.
#' @param K Number of realizations.
#'
#' @return Profile marginal log-likelihood.
marginal_ll_multi <- function(rho_idx, Y_mat, A, phi, n, K) {
  pr <- Q_precomp[[rho_idx]]

  Sy_shape <- diag(n) + phi * A %*% pr$Qinv %*% t(A)
  Sy_inv <- solve(Sy_shape)
  ld_Sy <- as.numeric(determinant(Sy_shape, logarithm = TRUE)$modulus)

  total_qform <- sum(diag(t(Y_mat) %*% Sy_inv %*% Y_mat))
  s2e <- total_qform / (n * K)

  if (!is.finite(s2e) || s2e <= 0) {
    return(-Inf)
  }

  -K * n / 2 * log(s2e) - K / 2 * ld_Sy - K * n / 2
}

# -------------------------------------------------------------------
# Marginal log-likelihood with phi profiled over a grid
# -------------------------------------------------------------------

#' Profile over phi for a fixed rho
#'
#' @param rho_idx Index into rho_grid / Q_precomp.
#' @param Y_mat n x K matrix of observations.
#' @param A Observation matrix.
#' @param n Number of observations per realization.
#' @param K Number of realizations.
#' @param phi_grid Candidate phi values.
#'
#' @return Best profile log-likelihood over phi_grid.
marginal_ll_profile_phi <- function(
    rho_idx,
    Y_mat,
    A,
    n,
    K,
    phi_grid = seq(0.5, 20, by = 0.5)
) {
  pr <- Q_precomp[[rho_idx]]
  AQinvAt <- A %*% pr$Qinv %*% t(A)

  best <- -Inf

  for (phi in phi_grid) {
    Sy_shape <- diag(n) + phi * AQinvAt

    Sy_inv <- tryCatch(
      solve(Sy_shape),
      error = function(e) NULL
    )

    if (is.null(Sy_inv)) {
      next
    }

    ld_Sy <- as.numeric(determinant(Sy_shape, logarithm = TRUE)$modulus)
    total_qf <- sum(diag(t(Y_mat) %*% Sy_inv %*% Y_mat))
    s2e <- total_qf / (n * K)

    if (!is.finite(s2e) || s2e <= 0) {
      next
    }

    ll <- -K * n / 2 * log(s2e) - K / 2 * ld_Sy - K * n / 2

    if (ll > best) {
      best <- ll
    }
  }

  best
}

# -------------------------------------------------------------------
# Plot output
# -------------------------------------------------------------------

pdf("vignettes/multi_realization_test.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# -------------------------------------------------------------------
# TEST 1: ll surface vs number of realizations
# -------------------------------------------------------------------

cat("--- Test 1: ll surface vs K (realizations) ---\n")

K_sizes <- c(1, 5, 20, 100, 500)
cols <- c("#2166ac", "#4dac26", "#e66101", "#762a83", "#d01c8b")
ll_all <- matrix(NA_real_, nrow = length(K_sizes), ncol = length(rho_grid))

for (ki in seq_along(K_sizes)) {
  K <- K_sizes[ki]
  set.seed(42)

  Y_mat <- simulate_Y_mat(
    K = K,
    Q_true = Q_true,
    A = A,
    phi_true = phi_true,
    sigma2e_true = sigma2e_true
  )

  ll_vals <- sapply(
    seq_along(rho_grid),
    marginal_ll_multi,
    Y_mat = Y_mat,
    A = A,
    phi = phi_true,
    n = n_per_field,
    K = K
  )

  ll_all[ki, ] <- ll_vals - max(ll_vals)
  peak <- rho_grid[which.max(ll_vals)]

  cat(sprintf("  K = %3d: peak rho = %.3f\n", K, peak))
}

plot(
  rho_grid, ll_all[1, ],
  type = "b", pch = 19, col = cols[1], cex = 0.6,
  xlab = "rho", ylab = "ll (normalized)",
  main = "Multi-realization marginal ll\nvs number of field realizations K",
  ylim = range(ll_all, na.rm = TRUE)
)

for (ki in 2:length(K_sizes)) {
  lines(
    rho_grid, ll_all[ki, ],
    type = "b", pch = 19, col = cols[ki], cex = 0.6
  )
}

abline(v = rho_true, col = "red", lty = 2, lwd = 2)
legend(
  "bottomleft",
  legend = sprintf("K=%d", K_sizes),
  col = cols,
  pch = 19,
  lty = 1,
  bty = "n",
  cex = 0.8
)

# -------------------------------------------------------------------
# TEST 2: replicated estimation across K realizations
# -------------------------------------------------------------------

cat("\n--- Test 2: replicated estimation ---\n")

n_reps <- 20L
results <- data.frame()

for (K in K_sizes) {
  rho_ests <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100)

    Y_mat <- simulate_Y_mat(
      K = K,
      Q_true = Q_true,
      A = A,
      phi_true = phi_true,
      sigma2e_true = sigma2e_true
    )

    ll_vals <- sapply(
      seq_along(rho_grid),
      marginal_ll_multi,
      Y_mat = Y_mat,
      A = A,
      phi = phi_true,
      n = n_per_field,
      K = K
    )

    rho_ests[rep] <- rho_grid[which.max(ll_vals)]
  }

  cat(sprintf(
    "  K = %3d: mean=%.3f sd=%.3f bias=%.3f\n",
    K, mean(rho_ests), sd(rho_ests), mean(rho_ests) - rho_true
  ))

  results <- rbind(
    results,
    data.frame(
      K = K,
      mean = mean(rho_ests),
      sd = sd(rho_ests),
      bias = mean(rho_ests) - rho_true,
      rmse = sqrt(mean((rho_ests - rho_true)^2))
    )
  )
}

plot(
  results$K, results$mean,
  type = "b", pch = 19, col = "#2166ac", log = "x",
  xlab = "K realizations (log)", ylab = "mean rho estimate",
  main = "Multi-realization: mean rho vs K",
  ylim = c(min(0.4, min(results$mean - results$sd)),
           max(1.0, max(results$mean + results$sd)))
)

abline(h = rho_true, col = "red", lty = 2, lwd = 2)
segments(
  results$K, results$mean - results$sd,
  results$K, results$mean + results$sd,
  col = "#2166ac", lwd = 2
)
legend(
  "topright",
  legend = sprintf("True=%.1f", rho_true),
  col = "red",
  lty = 2,
  bty = "n"
)

plot(
  results$K, results$sd,
  type = "b", pch = 19, col = "#e66101", log = "xy",
  xlab = "K realizations (log)", ylab = "SD (log)",
  main = "Multi-realization: SD vs K"
)

lines(
  results$K,
  results$sd[1] * sqrt(results$K[1] / results$K),
  col = "gray",
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c("Empirical SD", "1/sqrt(K)"),
  col = c("#e66101", "gray"),
  lty = c(1, 2),
  bty = "n"
)

# -------------------------------------------------------------------
# TEST 3: profile both phi and rho
# -------------------------------------------------------------------

cat("\n--- Test 3: profile both phi and rho ---\n")

K_sizes2 <- c(5, 20, 100)
cols2 <- c("#2166ac", "#4dac26", "#e66101")
ll_all2 <- matrix(NA_real_, nrow = length(K_sizes2), ncol = length(rho_grid))

for (ki in seq_along(K_sizes2)) {
  K <- K_sizes2[ki]
  set.seed(42)

  Y_mat <- simulate_Y_mat(
    K = K,
    Q_true = Q_true,
    A = A,
    phi_true = phi_true,
    sigma2e_true = sigma2e_true
  )

  ll_vals <- sapply(
    seq_along(rho_grid),
    marginal_ll_profile_phi,
    Y_mat = Y_mat,
    A = A,
    n = n_per_field,
    K = K
  )

  ll_all2[ki, ] <- ll_vals - max(ll_vals)
  peak <- rho_grid[which.max(ll_vals)]

  cat(sprintf("  K = %3d (phi profiled): peak rho = %.3f\n", K, peak))
}

plot(
  rho_grid, ll_all2[1, ],
  type = "b", pch = 19, col = cols2[1], cex = 0.6,
  xlab = "rho", ylab = "ll (normalized)",
  main = "Multi-realization, phi profiled\nvs K",
  ylim = range(ll_all2, na.rm = TRUE)
)

for (ki in 2:length(K_sizes2)) {
  lines(
    rho_grid, ll_all2[ki, ],
    type = "b", pch = 19, col = cols2[ki], cex = 0.6
  )
}

abline(v = rho_true, col = "red", lty = 2, lwd = 2)
legend(
  "bottomleft",
  legend = sprintf("K=%d", K_sizes2),
  col = cols2,
  pch = 19,
  lty = 1,
  bty = "n",
  cex = 0.8
)

# -------------------------------------------------------------------
# TEST 4: empirical vs theoretical covariance
# -------------------------------------------------------------------

cat("\n--- Test 4: empirical vs theoretical covariance ---\n")

set.seed(42)
K_large <- 200

Y_large <- simulate_Y_mat(
  K = K_large,
  Q_true = Q_true,
  A = A,
  phi_true = phi_true,
  sigma2e_true = sigma2e_true
)

# center rows before covariance estimation
row_means <- rowMeans(Y_large)
Y_centered <- Y_large - row_means
emp_cov <- tcrossprod(Y_centered) / K_large

idx_true <- which.min(abs(rho_grid - rho_true))
idx_09 <- which.min(abs(rho_grid - 0.9))

Q_t <- Q_precomp[[idx_true]]
Q_09 <- Q_precomp[[idx_09]]

theo_cov_true <- sigma2e_true * (
  diag(n_per_field) + phi_true * A %*% Q_t$Qinv %*% t(A)
)

Sy09 <- diag(n_per_field) + phi_true * A %*% Q_09$Qinv %*% t(A)
Sy09inv <- solve(Sy09)

s2e_09 <- sum(diag(t(Y_large) %*% Sy09inv %*% Y_large)) / (n_per_field * K_large)

theo_cov_09 <- s2e_09 * (
  diag(n_per_field) + phi_true * A %*% Q_09$Qinv %*% t(A)
)

plot(
  seq_len(n_per_field), emp_cov[1, ],
  type = "l", col = "black", lwd = 2,
  xlab = "observation j", ylab = "Cov(y_1, y_j)",
  main = "Test 4: empirical vs theoretical cov\n(row 1)"
)

lines(seq_len(n_per_field), theo_cov_true[1, ], col = "red", lwd = 2, lty = 2)
lines(seq_len(n_per_field), theo_cov_09[1, ], col = "#e66101", lwd = 2, lty = 3)

legend(
  "topright",
  legend = c(
    "Empirical",
    sprintf("Theo rho=%.1f", rho_grid[idx_true]),
    sprintf("Theo rho=%.1f (s2e=%.2f)", rho_grid[idx_09], s2e_09)
  ),
  col = c("black", "red", "#e66101"),
  lwd = 2,
  lty = c(1, 2, 3),
  bty = "n",
  cex = 0.8
)

par(mfrow = c(1, 1))
dev.off()

cat("Plots saved to vignettes/multi_realization_test.pdf\n")
