library(Matrix)

# ── Load project sources ─────────────────────────────────────────────────────
source("R/utils.R")
source("R/solvers.R")
source("R/fit.R")

# ── Build a SAR prior precision matrix on a grid ────────────────────────────
# Q_SAR = (I - rho * W)' (I - rho * W)  with row-standardised W
make_sar_Q <- function(nrow_grid, ncol_grid, rho = 0.95) {
  p   <- nrow_grid * ncol_grid
  # Queen adjacency on grid
  idx <- function(r, c) (c - 1) * nrow_grid + r
  i_idx <- integer(0); j_idx <- integer(0)
  for (r in seq_len(nrow_grid)) {
    for (c in seq_len(ncol_grid)) {
      nbrs <- list(c(r-1,c), c(r+1,c), c(r,c-1), c(r,c+1),
                   c(r-1,c-1), c(r-1,c+1), c(r+1,c-1), c(r+1,c+1))
      for (nb in nbrs) {
        if (nb[1] >= 1 && nb[1] <= nrow_grid && nb[2] >= 1 && nb[2] <= ncol_grid) {
          i_idx <- c(i_idx, idx(r, c))
          j_idx <- c(j_idx, idx(nb[1], nb[2]))
        }
      }
    }
  }
  # Row-standardised W
  counts <- tabulate(i_idx, nbins = p)
  w_vals <- 1 / counts[i_idx]
  W <- sparseMatrix(i = i_idx, j = j_idx, x = w_vals, dims = c(p, p))
  S <- Diagonal(p) - rho * W          # (I - rho W)
  Q <- Matrix::crossprod(S)            # S'S  -- symmetric positive definite
  Q
}

# ── Simulate data ────────────────────────────────────────────────────────────
set.seed(42)

GRID_ROWS <- 200    # p = 40x40 = 1600 grid cells
GRID_COLS <- 100
N_OBS     <- 5000   # n << p  (downscaling regime)
PHI       <- 10    # moderate-to-large signal-to-noise (hard case for PCG)
SIGMA2E   <- 1

cat(sprintf("Grid: %dx%d = %d cells | n_obs = %d | phi = %g\n\n",
            GRID_ROWS, GRID_COLS, GRID_ROWS * GRID_COLS, N_OBS, PHI))

p <- GRID_ROWS * GRID_COLS
Q <- make_sar_Q(GRID_ROWS, GRID_COLS, rho = 0.95)

# True field from SAR prior
R_chol <- chol(Q)
x_true <- as.numeric(solve(R_chol, rnorm(p)))

# Random sparse-ish observation matrix (each obs covers ~5 grid cells)
A <- matrix(0, N_OBS, p)
for (i in seq_len(N_OBS)) {
  cols <- sample(p, sample(3:7, 1))
  A[i, cols] <- runif(length(cols))
  A[i, ]     <- A[i, ] / sum(A[i, ])   # normalise rows
}
A <- Matrix(A, sparse = TRUE)

y <- as.numeric(A %*% x_true) + sqrt(SIGMA2E) * rnorm(N_OBS)

# ── Wire up operators ────────────────────────────────────────────────────────
apply_A    <- function(v) as.numeric(A %*% v)
apply_At   <- function(v) as.numeric(Matrix::crossprod(A, v))
apply_Q    <- function(v) as.numeric(Q %*% v)
apply_Rinv <- function(v) v    # R = I

apply_K <- make_apply_K(apply_A, apply_At, apply_Q, apply_Rinv, PHI)

# Preconditioner: sparse Cholesky of Q, apply Q^{-1}
CQ <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
apply_Qinv <- function(v) as.numeric(Matrix::solve(CQ, v))

# Preconditioned: M^{-1} = phi * Q^{-1}
precond_Q <- function(v) PHI * apply_Qinv(v)

# Right-hand side
rhs <- apply_At(apply_Rinv(y))

# ── Benchmark helper ─────────────────────────────────────────────────────────
run_bench <- function(label, precond, n_reps = 5) {
  iters <- numeric(n_reps)
  times <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    t0  <- proc.time()["elapsed"]
    res <- pcg(apply_K, rhs, precond = precond, tol = 1e-6, maxit = 10000L)
    times[r] <- proc.time()["elapsed"] - t0
    iters[r] <- res$iter
    if (!res$converged) warning(label, ": did not converge in rep ", r)
  }
  cat(sprintf("%-30s  iters: %4.0f (±%4.1f)   time: %.3fs (±%.4fs)  converged: %s\n",
              label,
              mean(iters), sd(iters),
              mean(times), sd(times),
              if (all(iters < 10000)) "YES" else "NO"))
  invisible(list(iters = iters, times = times))
}

cat("── PCG convergence benchmark ──────────────────────────────────────────\n")
cat(sprintf("  tol = 1e-6 | maxit = 10000 | %d reps each\n\n", 5L))

r_none  <- run_bench("No preconditioner",       precond = NULL)
r_precQ <- run_bench("Precond: phi * Q^{-1}",   precond = precond_Q)

cat("\n── Summary ────────────────────────────────────────────────────────────\n")
speedup_iter <- mean(r_none$iters)  / mean(r_precQ$iters)
speedup_time <- mean(r_none$times)  / mean(r_precQ$times)
cat(sprintf("  Iteration speedup:  %.1fx fewer iterations\n", speedup_iter))
cat(sprintf("  Wall-time speedup:  %.1fx faster\n", speedup_time))
cat(sprintf("  (Overhead of applying Q^{-1} per iter accounts for the difference)\n"))

# ── Sanity check: same solution? ─────────────────────────────────────────────
cat("\n── Solution agreement ─────────────────────────────────────────────────\n")
sol_none  <- pcg(apply_K, rhs, precond = NULL,      tol = 1e-8, maxit = 20000L)
sol_precQ <- pcg(apply_K, rhs, precond = precond_Q, tol = 1e-8, maxit = 20000L)
rel_err   <- sqrt(sum((sol_none$x - sol_precQ$x)^2)) / sqrt(sum(sol_none$x^2))
cat(sprintf("  Relative L2 diff between solutions: %.2e  %s\n",
            rel_err, if (rel_err < 1e-4) "(OK)" else "(MISMATCH - check!)"))

# ── Sensitivity to phi ───────────────────────────────────────────────────────
cat("\n── Iteration count vs phi (single run each) ───────────────────────────\n")
cat(sprintf("  %-8s  %10s  %10s  %8s\n", "phi", "iters_none", "iters_precQ", "speedup"))
make_precond <- function(phi) function(v) phi * apply_Qinv(v)

for (phi_test in c(0.1, 0.5, 1, 5, 10, 50, 100)) {
  K_test      <- make_apply_K(apply_A, apply_At, apply_Q, apply_Rinv, phi_test)
  pc_test     <- make_precond(phi_test)
  rhs_test    <- apply_At(y)   # same rhs, phi only enters K
  r_n <- pcg(K_test, rhs_test, precond = NULL,    tol = 1e-6, maxit = 20000L)
  r_p <- pcg(K_test, rhs_test, precond = pc_test, tol = 1e-6, maxit = 20000L)
  cat(sprintf("  %-8g  %10d  %10d  %8.1fx\n",
              phi_test, r_n$iter, r_p$iter, r_n$iter / max(r_p$iter, 1)))
}
