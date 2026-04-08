# Woodbury solver speed test
# Tests fit and posterior_se speed for p >> n regime
# Compares Woodbury vs Cholesky vs PCG where feasible
# Run with: pkgload::load_all("."); source("vignettes/woodbury_speed_test.R")

pkgload::load_all(".")
library(Matrix)

pdf("vignettes/woodbury_speed_test.pdf", width = 12, height = 10)

cat("=== Woodbury speed test (p >> n regime) ===\n\n")

# -----------------------------------------------------------------------
# Helper: make a p >> n problem
# Q = AR(1) precision so Q^{-1} is easy to apply via sparse solve
# -----------------------------------------------------------------------
make_problem_pn <- function(n, p, phi = 5, rho = 0.9, seed = 42) {
  set.seed(seed)
  A <- matrix(rnorm(n * p), n, p)
  Q <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                          diagonals = list(rep(-rho, p-1),
                                           c(1, rep(1 + rho^2, p-2), 1),
                                           rep(-rho, p-1)))
  Q <- Matrix::forceSymmetric(Q)

  # precompute sparse Cholesky of Q for fast Q^{-1} application
  CQ <- Matrix::Cholesky(Q)
  apply_Qinv <- function(v) as.numeric(Matrix::solve(CQ, v))

  x_true <- apply_Qinv(rnorm(p))
  y      <- as.numeric(A %*% x_true + rnorm(n))
  list(y = y, A = A, Q = Q, apply_Qinv = apply_Qinv, phi = phi, x_true = x_true, CQ = CQ)
}

# -----------------------------------------------------------------------
# TEST 1: Fit speed — Woodbury vs Cholesky vs PCG across p >> n sizes
# Fix n, vary p
# -----------------------------------------------------------------------
cat("--- Test 1: Fit speed, fixed n = 200, varying p ---\n")

n_fixed <- 200
p_sizes <- c(500, 1000, 2000, 5000, 10000)

results_fit <- data.frame()

for (p in p_sizes) {
  cat(sprintf("  p = %d\n", p))
  prob <- make_problem_pn(n_fixed, p)

  # Woodbury
  t_wb <- system.time(
    fit_wb <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi,
                          Q_inv = prob$apply_Qinv, solver = "woodbury")
  )["elapsed"]

  # PCG
  t_pcg <- system.time(
    fit_pcg <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
  )["elapsed"]

  # Cholesky only feasible for small p
  t_chol <- NA
  max_diff_chol <- NA
  if (p <= 2000) {
    t_chol <- system.time(
      fit_chol <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
    )["elapsed"]
    max_diff_chol <- max(abs(fit_wb$posterior_mean - fit_chol$posterior_mean))
  }

  max_diff_pcg <- max(abs(fit_wb$posterior_mean - fit_pcg$posterior_mean))

  cat(sprintf("    Woodbury: %.2fs  PCG: %.2fs  Cholesky: %s\n",
              t_wb, t_pcg, if (is.na(t_chol)) "skipped" else sprintf("%.2fs", t_chol)))
  cat(sprintf("    max_diff vs PCG: %.2e  vs Cholesky: %s\n",
              max_diff_pcg, if (is.na(max_diff_chol)) "skipped" else sprintf("%.2e", max_diff_chol)))

  results_fit <- rbind(results_fit, data.frame(
    p = p, n = n_fixed,
    t_wb = t_wb, t_pcg = t_pcg, t_chol = t_chol,
    max_diff_pcg = max_diff_pcg
  ))
}

# -----------------------------------------------------------------------
# TEST 2: Fit speed — fixed ratio p/n = 10, varying n
# -----------------------------------------------------------------------
cat("\n--- Test 2: Fit speed, fixed ratio p/n = 10 ---\n")

sizes_ratio <- list(
  c(n = 100,  p = 1000),
  c(n = 200,  p = 2000),
  c(n = 500,  p = 5000),
  c(n = 1000, p = 10000)
)

results_ratio <- data.frame()

for (sz in sizes_ratio) {
  n <- sz["n"]; p <- sz["p"]
  cat(sprintf("  n = %d, p = %d\n", n, p))
  prob <- make_problem_pn(n, p)

  t_wb <- system.time(
    fit_wb <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi,
                          Q_inv = prob$apply_Qinv, solver = "woodbury")
  )["elapsed"]

  t_pcg <- system.time(
    fit_pcg <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
  )["elapsed"]

  max_diff <- max(abs(fit_wb$posterior_mean - fit_pcg$posterior_mean))
  cat(sprintf("    Woodbury: %.2fs  PCG: %.2fs  max_diff: %.2e\n", t_wb, t_pcg, max_diff))

  results_ratio <- rbind(results_ratio, data.frame(
    n = n, p = p, t_wb = t_wb, t_pcg = t_pcg, max_diff = max_diff
  ))
}

# -----------------------------------------------------------------------
# TEST 3: posterior_se speed — Woodbury vs Cholesky
# -----------------------------------------------------------------------
cat("\n--- Test 3: posterior_se speed, fixed n = 200, varying p ---\n")

results_se <- data.frame()

for (p in c(500, 1000, 2000, 5000)) {
  cat(sprintf("  p = %d\n", p))
  prob    <- make_problem_pn(n_fixed, p)
  fit_wb  <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi,
                         Q_inv = prob$apply_Qinv, solver = "woodbury")

  t_se_wb <- system.time(se_wb <- posterior_se(fit_wb))["elapsed"]

  t_se_chol <- NA
  cor_se    <- NA
  if (p <= 2000) {
    fit_chol  <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
    t_se_chol <- system.time(se_chol <- posterior_se(fit_chol))["elapsed"]
    cor_se    <- cor(se_wb, se_chol)
    max_re    <- max(abs(se_wb - se_chol) / se_chol)
    cat(sprintf("    Woodbury SE: %.2fs  Cholesky SE: %.2fs  cor: %.4f  max_rel_err: %.2e\n",
                t_se_wb, t_se_chol, cor_se, max_re))
  } else {
    cat(sprintf("    Woodbury SE: %.2fs  Cholesky SE: skipped\n", t_se_wb))
  }

  results_se <- rbind(results_se, data.frame(
    p = p, t_se_wb = t_se_wb, t_se_chol = t_se_chol, cor_se = cor_se
  ))
}

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------

par(mfrow = c(2, 2))

# Plot 1: fit times, fixed n, varying p
plot(results_fit$p, results_fit$t_wb,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = sprintf("Fit time (n = %d, varying p)", n_fixed),
     ylim = range(c(results_fit$t_wb, results_fit$t_pcg), na.rm = TRUE))
lines(results_fit$p, results_fit$t_pcg, type = "b", pch = 19, col = "#e66101")
chol_rows <- !is.na(results_fit$t_chol)
if (any(chol_rows))
  lines(results_fit$p[chol_rows], results_fit$t_chol[chol_rows],
        type = "b", pch = 19, col = "#4dac26")
legend("topleft", legend = c("Woodbury", "PCG", "Cholesky"),
       col = c("#2166ac", "#e66101", "#4dac26"), pch = 19, lty = 1, bty = "n")

# Plot 2: fit times, fixed ratio p/n = 10
plot(results_ratio$p, results_ratio$t_wb,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "p (coefficients, p/n = 10)", ylab = "Time (s)",
     main = "Fit time (p/n = 10)",
     ylim = range(c(results_ratio$t_wb, results_ratio$t_pcg)))
lines(results_ratio$p, results_ratio$t_pcg, type = "b", pch = 19, col = "#e66101")
legend("topleft", legend = c("Woodbury", "PCG"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

# Plot 3: SE times
plot(results_se$p, results_se$t_se_wb,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = sprintf("posterior_se time (n = %d)", n_fixed),
     ylim = range(c(results_se$t_se_wb, results_se$t_se_chol), na.rm = TRUE))
chol_se_rows <- !is.na(results_se$t_se_chol)
if (any(chol_se_rows))
  lines(results_se$p[chol_se_rows], results_se$t_se_chol[chol_se_rows],
        type = "b", pch = 19, col = "#4dac26")
legend("topleft", legend = c("Woodbury SE", "Cholesky SE"),
       col = c("#2166ac", "#4dac26"), pch = 19, lty = 1, bty = "n")

# Plot 4: agreement Woodbury vs PCG
plot(results_fit$p, results_fit$max_diff_pcg,
     type = "b", pch = 19, col = "#762a83", log = "y",
     xlab = "p (coefficients)", ylab = "Max |Woodbury - PCG|",
     main = "Woodbury vs PCG agreement")

par(mfrow = c(1, 1))

cat("\n=== Woodbury speed test complete ===\n")
dev.off()
cat("Plots saved to vignettes/woodbury_speed_test.pdf\n")
