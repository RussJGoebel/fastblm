# =============================================================================
# test_rho_recovery.R
#
# Model assumed by tune_reml():
#   y   = A x + eps,   eps ~ N(0, sigma2e * I)
#   x   ~ N(0, phi * sigma2e * Q^{-1})
#
# where Q = (I - rho*W)'(I - rho*W),  W row-normalised adjacency.
#
# To simulate from this:
#   x = sqrt(phi * sigma2e) * (I - rho*W)^{-1} z,   z ~ N(0, I)
#   y = A x + sqrt(sigma2e) * noise,                 noise ~ N(0, I)
#
# Scenarios:
#   A) A = I  (direct, one obs per pixel)
#   B) A = footprint average matrix, moderate aggregation
#   C) A = footprint average matrix, heavy aggregation
#
# Usage:
#   source("utils.R"); source("solvers.R"); source("fit.R"); source("reml.R")
#   source("test_rho_recovery.R")
# =============================================================================

set.seed(42)

# ---- simulation helpers -----------------------------------------------------

#' Row-normalised queen-adjacency W on an m x m grid
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

#' Simulate x from SAR prior:  x ~ N(0, phi * sigma2e * Q^{-1})
#'   Q = (I - rho*W)'(I - rho*W)
#'   => x = sqrt(phi * sigma2e) * (I - rho*W)^{-1} z,  z ~ N(0, I)
simulate_sar <- function(W, rho, phi, sigma2e) {
  p <- nrow(W)
  M <- Matrix::Diagonal(p) - rho * W          # (I - rho*W)
  z <- rnorm(p)
  as.numeric(Matrix::solve(M, z)) * sqrt(phi * sigma2e)
}

#' Sliding-window footprint-average matrix (n_obs x p)
make_agg_matrix <- function(p, n_obs, fsize) {
  stride <- max(1L, floor(p / n_obs))
  starts <- seq(1L, by = stride, length.out = n_obs)
  ends   <- pmin(starts + fsize - 1L, p)
  rows   <- rep(seq_len(n_obs), times = ends - starts + 1L)
  cols   <- unlist(mapply(seq, starts, ends, SIMPLIFY = FALSE))
  vals   <- rep(1 / (ends - starts + 1L), times = ends - starts + 1L)
  Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_obs, p))
}

# ---- REML helper ------------------------------------------------------------

#' Q_fun factory for tune_reml
#' theta = c(rho = ...);  Q = (I - rho*W)'(I - rho*W)
make_sar_Q_fun <- function(W) {
  p <- nrow(W)
  function(theta) {
    rho <- theta[["rho"]]
    M   <- Matrix::Diagonal(p) - rho * W
    Q   <- Matrix::crossprod(M)
    # log|Q| = log|(I-rhoW)'(I-rhoW)| = 2 * log|det(I-rhoW)|
    ld  <- 2 * as.numeric(Matrix::determinant(M, logarithm = TRUE)$modulus)
    list(Q = Q, log_det_Q = ld, Q_matrix = Q)
  }
}

#' Run one recovery experiment
run_recovery <- function(label, y, A, W, true_rho, true_phi, verbose = FALSE) {
  Q_fun <- make_sar_Q_fun(W)
  res   <- tune_reml(
    y             = y,
    A             = as.matrix(A),
    Q_fun         = Q_fun,
    theta_init    = c(rho = 0.5),
    lower         = c(rho = 0.0),
    upper         = c(rho = 0.999),
    logdet_method = "cholesky",
    verbose       = verbose
  )
  data.frame(
    scenario  = label,
    true_rho  = true_rho,
    true_phi  = true_phi,
    est_rho   = round(res$theta[["rho"]], 4),
    est_phi   = round(res$phi, 4),
    rho_error = round(abs(res$theta[["rho"]] - true_rho), 4),
    converged = (res$optim$convergence == 0L),
    ll        = round(res$value, 2)
  )
}

# =============================================================================
# Setup
# =============================================================================
cat("=====================================================\n")
cat(" SAR rho recovery test\n")
cat("=====================================================\n\n")

m       <- 20        # 20x20 = 400 latent pixels
p       <- m * m
phi     <- 5.0       # true signal-to-noise ratio (sigma2b / sigma2e)
sigma2e <- 1.0       # true noise variance

cat(sprintf("Grid: %dx%d = %d pixels\n", m, m, p))
cat(sprintf("True phi = %.1f,  sigma2e = %.1f\n", phi, sigma2e))
cat(sprintf("=> sigma2b = phi * sigma2e = %.1f\n\n", phi * sigma2e))

W       <- make_grid_W(m)
results <- list()

# =============================================================================
# Scenario A: direct (A = I, n = p)
# =============================================================================
cat("--- Scenario A: direct observations (A = I, n = p = 400) ---\n\n")

for (true_rho in c(0.3, 0.6, 0.9)) {
  cat(sprintf("  rho = %.1f ...\n", true_rho))
  x   <- simulate_sar(W, true_rho, phi, sigma2e)
  y   <- x + rnorm(p, sd = sqrt(sigma2e))

  res <- run_recovery(
    label    = sprintf("A_direct   rho=%.1f", true_rho),
    y        = y,
    A        = Matrix::Diagonal(p),
    W        = W,
    true_rho = true_rho,
    true_phi = phi
  )
  print(res, row.names = FALSE); cat("\n")
  results[[length(results) + 1]] <- res
}

# =============================================================================
# Scenario B: moderate aggregation (n=150, footprint=8px)
# =============================================================================
cat("--- Scenario B: moderate aggregation (n=150, footprint=8px) ---\n\n")

A_B <- make_agg_matrix(p, n_obs = 150L, fsize = 8L)

for (true_rho in c(0.3, 0.6, 0.9)) {
  cat(sprintf("  rho = %.1f ...\n", true_rho))
  x   <- simulate_sar(W, true_rho, phi, sigma2e)
  y   <- as.numeric(A_B %*% x) + rnorm(nrow(A_B), sd = sqrt(sigma2e))

  res <- run_recovery(
    label    = sprintf("B_moderate rho=%.1f", true_rho),
    y        = y,
    A        = A_B,
    W        = W,
    true_rho = true_rho,
    true_phi = phi
  )
  print(res, row.names = FALSE); cat("\n")
  results[[length(results) + 1]] <- res
}

# =============================================================================
# Scenario C: heavy aggregation (n=50, footprint=30px)
# =============================================================================
cat("--- Scenario C: heavy aggregation (n=50, footprint=30px) ---\n\n")

A_C <- make_agg_matrix(p, n_obs = 50L, fsize = 30L)

for (true_rho in c(0.3, 0.9)) {
  cat(sprintf("  rho = %.1f ...\n", true_rho))
  x   <- simulate_sar(W, true_rho, phi, sigma2e)
  y   <- as.numeric(A_C %*% x) + rnorm(nrow(A_C), sd = sqrt(sigma2e))

  res <- run_recovery(
    label    = sprintf("C_heavy    rho=%.1f", true_rho),
    y        = y,
    A        = A_C,
    W        = W,
    true_rho = true_rho,
    true_phi = phi
  )
  print(res, row.names = FALSE); cat("\n")
  results[[length(results) + 1]] <- res
}

# =============================================================================
# Summary
# =============================================================================
cat("=====================================================\n")
cat(" SUMMARY\n")
cat("=====================================================\n")
print(do.call(rbind, results), row.names = FALSE)

cat("\nDiagnosis:\n")
cat("  rho_error < 0.10  ->  good recovery\n")
cat("  rho_error > 0.20  ->  poor recovery\n")
cat("  A recovers, B/C don't  ->  aggregation kills identifiability (expected)\n")
cat("  A also fails           ->  likely a code issue worth investigating\n")
