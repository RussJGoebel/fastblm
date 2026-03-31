# Sparse vs dense posterior_se speed test
# Tests whether sparse A dramatically speeds up posterior_se
# Run with: pkgload::load_all("."); source("vignettes/sparse_se_test.R")

pkgload::load_all(".")
library(Matrix)

cat("=== Sparse vs dense posterior_se speed test ===\n\n")

# -----------------------------------------------------------------------
# Helper: make sparse A (e.g. each row has only k nonzeros)
# This mimics a spatial integration matrix where each observation
# only touches a few basis functions
# -----------------------------------------------------------------------
make_sparse_A <- function(n, p, nnz_per_row = 5, seed = 42) {
  set.seed(seed)
  # each row has nnz_per_row nonzero entries at random columns
  rows <- rep(seq_len(n), each = nnz_per_row)
  cols <- unlist(lapply(seq_len(n), function(i) sample(p, nnz_per_row)))
  vals <- rnorm(n * nnz_per_row)
  Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, p))
}

make_dense_A <- function(n, p, seed = 42) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

make_Q <- function(p, rho = 0.5) {
  Q <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                          diagonals = list(rep(-rho, p-1),
                                           c(1, rep(1 + rho^2, p-2), 1),
                                           rep(-rho, p-1)))
  Matrix::forceSymmetric(Q)
}

# -----------------------------------------------------------------------
# TEST: sparse vs dense A across problem sizes
# -----------------------------------------------------------------------
cat(sprintf("%-6s %-6s | %-12s %-12s | %-12s %-12s | %s\n",
            "p", "n", "dense fit", "dense SE", "sparse fit", "sparse SE", "SE speedup"))
cat(strrep("-", 80), "\n")

p_sizes <- c(500, 1000, 2000, 5000)
phi <- 5

results <- data.frame()

for (p in p_sizes) {
  n <- max(50, p %/% 5)
  Q <- make_Q(p)
  y <- rnorm(n)

  # --- dense A ---
  A_dense <- make_dense_A(n, p)
  t_fit_dense <- system.time(
    fit_dense <- fit_fastblm(y, A_dense, Q, phi = phi, solver = "cholesky")
  )["elapsed"]
  t_se_dense <- system.time(
    se_dense <- posterior_se(fit_dense)
  )["elapsed"]

  # --- sparse A (5 nonzeros per row) ---
  A_sparse <- make_sparse_A(n, p, nnz_per_row = 5)
  t_fit_sparse <- system.time(
    fit_sparse <- fit_fastblm(y, A_sparse, Q, phi = phi, solver = "cholesky")
  )["elapsed"]
  t_se_sparse <- system.time(
    se_sparse <- posterior_se(fit_sparse)
  )["elapsed"]

  # --- constrain with q=1, sparse A ---
  C <- matrix(rnorm(p), nrow = 1)
  t_con_sparse <- system.time({
    fit_sparse_c <- constrain(fit_sparse, C)
    posterior_se(fit_sparse_c)
  })["elapsed"]

  # --- constrain with q=1, dense A ---
  t_con_dense <- system.time({
    fit_dense_c <- constrain(fit_dense, C)
    posterior_se(fit_dense_c)
  })["elapsed"]

  speedup    <- t_se_dense / t_se_sparse
  speedup_con <- t_con_dense / t_con_sparse

  cat(sprintf("%-6d %-6d | %-12.3f %-12.3f | %-12.3f %-12.3f | %.1fx  (constrain: dense=%.3fs sparse=%.3fs speedup=%.1fx)\n",
              p, n, t_fit_dense, t_se_dense, t_fit_sparse, t_se_sparse, speedup,
              t_con_dense, t_con_sparse, speedup_con))

  results <- rbind(results, data.frame(
    p = p, n = n,
    t_fit_dense = t_fit_dense, t_se_dense = t_se_dense,
    t_fit_sparse = t_fit_sparse, t_se_sparse = t_se_sparse,
    t_con_dense = t_con_dense, t_con_sparse = t_con_sparse,
    speedup = speedup, speedup_con = speedup_con
  ))
}

# -----------------------------------------------------------------------
# TEST 2: vary sparsity of A (nnz per row), fixed p
# -----------------------------------------------------------------------
cat("\n--- Varying sparsity of A, p = 2000 ---\n")
p <- 2000
n <- 400
Q <- make_Q(p)
y <- rnorm(n)

cat(sprintf("%-15s | %-12s %-12s\n", "nnz_per_row", "fit time", "SE time"))
cat(strrep("-", 45), "\n")

results2 <- data.frame()
for (nnz in c(2, 5, 10, 20, 50, 100, p)) {
  label <- if (nnz == p) "dense" else as.character(nnz)
  A <- if (nnz == p) make_dense_A(n, p) else make_sparse_A(n, p, nnz_per_row = nnz)
  t_fit <- system.time(fit <- fit_fastblm(y, A, Q, phi = phi, solver = "cholesky"))["elapsed"]
  t_se  <- system.time(posterior_se(fit))["elapsed"]
  cat(sprintf("%-15s | %-12.3f %-12.3f\n", label, t_fit, t_se))
  results2 <- rbind(results2, data.frame(nnz=nnz, t_fit=t_fit, t_se=t_se))
}

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------
pdf("vignettes/sparse_se_test.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

# Plot 1: SE time sparse vs dense
plot(results$p, results$t_se_dense,
     type = "b", pch = 19, col = "#e66101",
     xlab = "p", ylab = "posterior_se time (s)",
     main = "posterior_se: sparse vs dense A",
     ylim = range(c(results$t_se_dense, results$t_se_sparse)))
lines(results$p, results$t_se_sparse, type = "b", pch = 19, col = "#2166ac")
legend("topleft", legend = c("dense A", "sparse A (5 nnz/row)"),
       col = c("#e66101", "#2166ac"), pch = 19, lty = 1, bty = "n")

# Plot 2: SE speedup
plot(results$p, results$speedup,
     type = "b", pch = 19, col = "#4dac26",
     xlab = "p", ylab = "Speedup (dense / sparse)",
     main = "SE speedup from sparse A")
abline(h = 1, col = "red", lty = 2)

# Plot 3: constrain + SE time sparse vs dense
plot(results$p, results$t_con_dense,
     type = "b", pch = 19, col = "#e66101",
     xlab = "p", ylab = "constrain + posterior_se time (s)",
     main = "constrain + posterior_se: sparse vs dense A\n(q = 1)",
     ylim = range(c(results$t_con_dense, results$t_con_sparse)))
lines(results$p, results$t_con_sparse, type = "b", pch = 19, col = "#2166ac")
legend("topleft", legend = c("dense A", "sparse A (5 nnz/row)"),
       col = c("#e66101", "#2166ac"), pch = 19, lty = 1, bty = "n")

# Plot 4: constrain speedup
plot(results$p, results$speedup_con,
     type = "b", pch = 19, col = "#762a83",
     xlab = "p", ylab = "Speedup (dense / sparse)",
     main = "constrain + SE speedup from sparse A")
abline(h = 1, col = "red", lty = 2)

par(mfrow = c(1, 1))
dev.off()

cat("\nPlots saved to vignettes/sparse_se_test.pdf\n")
cat("=== Done ===\n")
