# =============================================================================
# test_rho_recovery.R
#
# Model assumed by tune_reml():
#   y   = A x + eps,   eps ~ N(0, sigma2e * I)
#   x   ~ N(0, phi * sigma2e * Q^{-1})
#   Q   = (I - rho*W)'(I - rho*W),  W row-normalised queen adjacency
#
# Simulation:
#   x = sqrt(phi * sigma2e) * (I - rho*W)^{-1} z,  z ~ N(0, I)
#   y = A x + sqrt(sigma2e) * noise
#
# Scenarios:
#   A) Multiple direct obs per pixel (n = k*p >> p). Well-posed baseline.
#      Each pixel observed k=5 times with independent noise.
#   B) Footprint averages, moderate aggregation (n=150, fsize=8)
#   C) Footprint averages, heavy aggregation  (n=50,  fsize=30)
#
# Scenario A tests that the code works. If A fails, there's a bug.
# Scenarios B/C test identifiability under aggregation.
#
# Usage:
#   source("utils.R"); source("solvers.R"); source("fit.R"); source("reml.R")
#   source("test_rho_recovery.R")
# =============================================================================

set.seed(42)

# ---- helpers ----------------------------------------------------------------

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

#' Simulate x ~ N(0, phi * sigma2e * Q^{-1}),  Q = (I-rhoW)'(I-rhoW)
simulate_sar <- function(W, rho, phi, sigma2e) {
  p <- nrow(W)
  M <- Matrix::Diagonal(p) - rho * W
  as.numeric(Matrix::solve(M, rnorm(p))) * sqrt(phi * sigma2e)
}

#' Replicated direct design: stack k copies of I_p  (n = k*p, A = rep(I, k))
make_replicated_A <- function(p, k) {
  # k independent observations per pixel => A is k*p x p block of identities
  do.call(rbind, replicate(k, diag(p), simplify = FALSE))
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

#' Q_fun factory: theta = c(rho = ...),  Q = (I-rhoW)'(I-rhoW)
make_sar_Q_fun <- function(W) {
  p <- nrow(W)
  function(theta) {
    rho <- theta[["rho"]]
    M   <- Matrix::Diagonal(p) - rho * W
    Q   <- Matrix::crossprod(M)
    ld  <- 2 * as.numeric(Matrix::determinant(M, logarithm = TRUE)$modulus)
    list(Q = Q, log_det_Q = ld, Q_matrix = Q)
  }
}

#' Run one recovery experiment, return tidy one-row data.frame
run_recovery <- function(label, y, A, W, true_rho, true_phi, verbose = FALSE) {
  res <- tune_reml(
    y             = y,
    A             = as.matrix(A),
    Q_fun         = make_sar_Q_fun(W),
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
phi     <- 5.0
sigma2e <- 1.0

W       <- make_grid_W(m)
results <- list()

cat(sprintf("Grid: %dx%d = %d pixels,  phi=%.1f,  sigma2e=%.1f\n\n",
            m, m, p, phi, sigma2e))

# =============================================================================
# Scenario A: replicated direct observations (n = 5*p = 2000 >> p = 400)
# Each pixel observed 5 times independently. n > p so phi is identifiable.
# This is the well-posed baseline: if rho fails here, it's a bug.
# =============================================================================
cat("--- Scenario A: replicated direct obs (k=5 reps/pixel, n=2000) ---\n\n")

k   <- 5L
A_A <- make_replicated_A(p, k)   # 2000 x 400

for (true_rho in c(0.3, 0.6, 0.9)) {
  cat(sprintf("  rho = %.1f ...\n", true_rho))
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_A %*% x) + rnorm(nrow(A_A), sd = sqrt(sigma2e))

  res <- run_recovery(
    label    = sprintf("A_direct   rho=%.1f", true_rho),
    y        = y,
    A        = A_A,
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
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_B %*% x) + rnorm(nrow(A_B), sd = sqrt(sigma2e))

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
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_C %*% x) + rnorm(nrow(A_C), sd = sqrt(sigma2e))

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
cat("  A fails (rho_error > 0.10)  ->  code bug\n")
cat("  A ok, B/C degrade           ->  aggregation kills identifiability (expected)\n")
cat("  rho=0.9 recovers better than rho=0.3 in B/C  ->  also expected\n")# =============================================================================
# test_rho_recovery.R
#
# Model assumed by tune_reml():
#   y   = A x + eps,   eps ~ N(0, sigma2e * I)
#   x   ~ N(0, phi * sigma2e * Q^{-1})
#   Q   = (I - rho*W)'(I - rho*W),  W row-normalised queen adjacency
#
# Simulation:
#   x = sqrt(phi * sigma2e) * (I - rho*W)^{-1} z,  z ~ N(0, I)
#   y = A x + sqrt(sigma2e) * noise
#
# Scenarios:
#   A) Multiple direct obs per pixel (n = k*p >> p). Well-posed baseline.
#      Each pixel observed k=5 times with independent noise.
#   B) Footprint averages, moderate aggregation (n=150, fsize=8)
#   C) Footprint averages, heavy aggregation  (n=50,  fsize=30)
#
# Scenario A tests that the code works. If A fails, there's a bug.
# Scenarios B/C test identifiability under aggregation.
#
# Usage:
#   source("utils.R"); source("solvers.R"); source("fit.R"); source("reml.R")
#   source("test_rho_recovery.R")
# =============================================================================

set.seed(42)

# ---- helpers ----------------------------------------------------------------

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

#' Simulate x ~ N(0, phi * sigma2e * Q^{-1}),  Q = (I-rhoW)'(I-rhoW)
simulate_sar <- function(W, rho, phi, sigma2e) {
  p <- nrow(W)
  M <- Matrix::Diagonal(p) - rho * W
  as.numeric(Matrix::solve(M, rnorm(p))) * sqrt(phi * sigma2e)
}

#' Replicated direct design: stack k copies of I_p  (n = k*p, A = rep(I, k))
make_replicated_A <- function(p, k) {
  # k independent observations per pixel => A is k*p x p block of identities
  do.call(rbind, replicate(k, diag(p), simplify = FALSE))
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

#' Q_fun factory: theta = c(rho = ...),  Q = (I-rhoW)'(I-rhoW)
make_sar_Q_fun <- function(W) {
  p <- nrow(W)
  function(theta) {
    rho <- theta[["rho"]]
    M   <- Matrix::Diagonal(p) - rho * W
    Q   <- Matrix::crossprod(M)
    ld  <- 2 * as.numeric(Matrix::determinant(M, logarithm = TRUE)$modulus)
    list(Q = Q, log_det_Q = ld, Q_matrix = Q)
  }
}

#' Run one recovery experiment, return tidy one-row data.frame
run_recovery <- function(label, y, A, W, true_rho, true_phi, verbose = FALSE) {
  res <- tune_reml(
    y             = y,
    A             = as.matrix(A),
    Q_fun         = make_sar_Q_fun(W),
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
phi     <- 5.0
sigma2e <- 1.0

W       <- make_grid_W(m)
results <- list()

cat(sprintf("Grid: %dx%d = %d pixels,  phi=%.1f,  sigma2e=%.1f\n\n",
            m, m, p, phi, sigma2e))

# =============================================================================
# Scenario A: replicated direct observations (n = 5*p = 2000 >> p = 400)
# Each pixel observed 5 times independently. n > p so phi is identifiable.
# This is the well-posed baseline: if rho fails here, it's a bug.
# =============================================================================
cat("--- Scenario A: replicated direct obs (k=5 reps/pixel, n=2000) ---\n\n")

k   <- 5L
A_A <- make_replicated_A(p, k)   # 2000 x 400

for (true_rho in c(0.3, 0.6, 0.9)) {
  cat(sprintf("  rho = %.1f ...\n", true_rho))
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_A %*% x) + rnorm(nrow(A_A), sd = sqrt(sigma2e))

  res <- run_recovery(
    label    = sprintf("A_direct   rho=%.1f", true_rho),
    y        = y,
    A        = A_A,
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
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_B %*% x) + rnorm(nrow(A_B), sd = sqrt(sigma2e))

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
  x <- simulate_sar(W, true_rho, phi, sigma2e)
  y <- as.numeric(A_C %*% x) + rnorm(nrow(A_C), sd = sqrt(sigma2e))

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
cat("  A fails (rho_error > 0.10)  ->  code bug\n")
cat("  A ok, B/C degrade           ->  aggregation kills identifiability (expected)\n")
cat("  rho=0.9 recovers better than rho=0.3 in B/C  ->  also expected\n")
