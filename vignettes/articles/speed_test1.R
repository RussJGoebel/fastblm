# fastblm speed tests
# Goal: ~5 minutes total runtime
# Tests:
#   1. Cholesky vs PCG fit across problem sizes
#   2. posterior_se: Cholesky vs Hutchinson speed and accuracy
#   3. Hutchinson convergence as n_probes increases

library(Matrix)
pkgload::load_all(".")

pdf("vignettes/speed_tests.pdf", width = 12, height = 10)

cat("=== fastblm speed tests ===\n\n")

# -----------------------------------------------------------------------
# Helper: make a test problem of size (n, p) with sparse AR(1) Q
# -----------------------------------------------------------------------
make_problem <- function(n, p, phi = 5, rho = 0.9, seed = 42) {
  set.seed(seed)
  A       <- matrix(rnorm(n * p), n, p)
  Q       <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                                diagonals = list(rep(-rho, p-1),
                                                 c(1, rep(1 + rho^2, p-2), 1),
                                                 rep(-rho, p-1)))
  Q       <- Matrix::forceSymmetric(Q)
  x_true  <- as.numeric(Matrix::solve(Matrix::chol(Q), rnorm(p)))
  y       <- as.numeric(A %*% x_true + rnorm(n))
  list(y = y, A = A, Q = Q, phi = phi, x_true = x_true)
}

# -----------------------------------------------------------------------
# TEST 1: Cholesky vs PCG fit across problem sizes
# Cholesky is O(p^3) so capped at p=2000
# PCG scales much better so goes further
# -----------------------------------------------------------------------
cat("--- Test 1: Cholesky vs PCG fit times ---\n")

sizes_chol <- list(
  c(n = 500,  p = 100),
  c(n = 1000, p = 500),
  c(n = 2000, p = 1000),
  c(n = 5000, p = 2000)
)

sizes_pcg <- list(
  c(n = 500,   p = 100),
  c(n = 1000,  p = 500),
  c(n = 2000,  p = 1000),
  c(n = 5000,  p = 2000),
  c(n = 10000, p = 5000),
  c(n = 50000, p = 10000)
)

# --- Cholesky fits ---
cat("  Cholesky:\n")
results_chol <- data.frame()
for (sz in sizes_chol) {
  n <- sz["n"]; p <- sz["p"]
  prob  <- make_problem(n, p)
  t     <- system.time(
    fit_c <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
  )["elapsed"]
  cor_true <- cor(fit_c$posterior_mean, prob$x_true)
  cat(sprintf("    n = %5d, p = %5d: %.2fs  cor(x_true): %.3f\n", n, p, t, cor_true))
  results_chol <- rbind(results_chol, data.frame(n=n, p=p, time=t, cor_true=cor_true))
}

# --- PCG fits ---
cat("  PCG:\n")
results_pcg <- data.frame()
for (sz in sizes_pcg) {
  n <- sz["n"]; p <- sz["p"]
  prob  <- make_problem(n, p)
  t     <- system.time(
    fit_p <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
  )["elapsed"]
  cor_true <- cor(fit_p$posterior_mean, prob$x_true)
  cat(sprintf("    n = %5d, p = %5d: %.2fs  cor(x_true): %.3f\n", n, p, t, cor_true))
  results_pcg <- rbind(results_pcg, data.frame(n=n, p=p, time=t, cor_true=cor_true))
}

# --- Agreement check on shared sizes ---
cat("  Agreement (Cholesky vs PCG, shared sizes):\n")
for (sz in sizes_chol) {
  n <- sz["n"]; p <- sz["p"]
  prob  <- make_problem(n, p)
  fit_c <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
  fit_p <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
  cat(sprintf("    n = %5d, p = %5d: max_diff = %.2e\n", n, p,
              max(abs(fit_c$posterior_mean - fit_p$posterior_mean))))
}

# Plot 1: fit times
par(mfrow = c(1, 2))
plot(results_pcg$p, results_pcg$time,
     type = "b", pch = 19, col = "#e66101",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = "Fit time: Cholesky vs PCG",
     ylim = range(c(results_chol$time, results_pcg$time)))
lines(results_chol$p, results_chol$time, type = "b", pch = 19, col = "#2166ac")
legend("topleft", legend = c("Cholesky", "PCG"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

plot(results_pcg$p, results_pcg$cor_true,
     type = "b", pch = 19, col = "#e66101",
     xlab = "p (coefficients)", ylab = "cor(posterior mean, x_true)",
     main = "Recovery of true x",
     ylim = c(0, 1))
lines(results_chol$p, results_chol$cor_true, type = "b", pch = 19, col = "#2166ac")
legend("bottomleft", legend = c("Cholesky", "PCG"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------
# TEST 2: posterior_se speed — Cholesky vs Hutchinson, across sizes
# -----------------------------------------------------------------------
cat("\n--- Test 2: posterior_se Cholesky vs Hutchinson speed ---\n")

# use Cholesky-feasible sizes for exact SE comparison
se_sizes <- list(
  c(n = 500,  p = 100),
  c(n = 1000, p = 500),
  c(n = 2000, p = 1000),
  c(n = 5000, p = 2000)
)

results_se <- data.frame()

for (sz in se_sizes) {
  n <- sz["n"]; p <- sz["p"]
  cat(sprintf("  n = %d, p = %d\n", n, p))
  prob  <- make_problem(n, p)
  fit_c <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")

  # Cholesky SE (exact)
  t_se_chol <- system.time({
    se_chol <- posterior_se(fit_c)
  })["elapsed"]

  # Hutchinson SE via PCG fit
  fit_p <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
  t_se_hutch <- system.time({
    se_hutch <- posterior_se(fit_p, n_probes = 100L)
  })["elapsed"]

  cor_se  <- cor(se_chol, se_hutch)
  rel_err <- mean(abs(se_hutch - se_chol) / se_chol)

  cat(sprintf("    Cholesky SE: %.2fs  Hutchinson SE: %.2fs  cor: %.3f  mean rel err: %.3f\n",
              t_se_chol, t_se_hutch, cor_se, rel_err))

  results_se <- rbind(results_se, data.frame(
    n = n, p = p,
    t_se_chol = t_se_chol, t_se_hutch = t_se_hutch,
    cor_se = cor_se, rel_err = rel_err
  ))
}

# Plot 2: SE times and accuracy
par(mfrow = c(1, 2))
plot(results_se$p, results_se$t_se_chol,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = "posterior_se time: Cholesky vs Hutchinson",
     ylim = range(c(results_se$t_se_chol, results_se$t_se_hutch)))
lines(results_se$p, results_se$t_se_hutch, type = "b", pch = 19, col = "#e66101")
legend("topleft", legend = c("Cholesky", "Hutchinson (100 probes)"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

plot(results_se$p, results_se$cor_se,
     type = "b", pch = 19, col = "#4dac26",
     xlab = "p (coefficients)", ylab = "Correlation",
     main = "Hutchinson vs Cholesky SE correlation",
     ylim = c(0, 1))
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------
# TEST 3: Hutchinson convergence as n_probes increases
# -----------------------------------------------------------------------
cat("\n--- Test 3: Hutchinson convergence with n_probes ---\n")

# fix a moderate size problem
n <- 1000; p <- 500
prob  <- make_problem(n, p)
fit_c <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
fit_p <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "pcg")
se_exact <- posterior_se(fit_c)

probe_counts <- c(10, 25, 50, 100, 200, 500)
results_probes <- data.frame()

for (np in probe_counts) {
  cat(sprintf("  n_probes = %d\n", np))
  # run 5 replicates to assess variance
  cors <- rel_errs <- times <- numeric(5)
  for (rep in 1:5) {
    t <- system.time({
      se_h <- posterior_se(fit_p, n_probes = np)
    })["elapsed"]
    cors[rep]     <- cor(se_h, se_exact)
    rel_errs[rep] <- mean(abs(se_h - se_exact) / se_exact)
    times[rep]    <- t
  }
  cat(sprintf("    cor: %.3f (sd %.4f)  rel_err: %.3f (sd %.4f)  time: %.2fs\n",
              mean(cors), sd(cors), mean(rel_errs), sd(rel_errs), mean(times)))

  results_probes <- rbind(results_probes, data.frame(
    n_probes  = np,
    cor_mean  = mean(cors),
    cor_sd    = sd(cors),
    err_mean  = mean(rel_errs),
    err_sd    = sd(rel_errs),
    time_mean = mean(times)
  ))
}

# Plot 3: Hutchinson convergence
par(mfrow = c(1, 3))

# correlation vs n_probes
plot(results_probes$n_probes, results_probes$cor_mean,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "n_probes", ylab = "Correlation with exact SE",
     main = "Hutchinson convergence: correlation",
     ylim = c(0, 1))
segments(results_probes$n_probes,
         results_probes$cor_mean - results_probes$cor_sd,
         results_probes$n_probes,
         results_probes$cor_mean + results_probes$cor_sd,
         col = "#2166ac")

# relative error vs n_probes
plot(results_probes$n_probes, results_probes$err_mean,
     type = "b", pch = 19, col = "#e66101",
     xlab = "n_probes", ylab = "Mean relative error",
     main = "Hutchinson convergence: relative error")
segments(results_probes$n_probes,
         results_probes$err_mean - results_probes$err_sd,
         results_probes$n_probes,
         results_probes$err_mean + results_probes$err_sd,
         col = "#e66101")

# time vs n_probes
plot(results_probes$n_probes, results_probes$time_mean,
     type = "b", pch = 19, col = "#4dac26",
     xlab = "n_probes", ylab = "Time (s)",
     main = "Hutchinson time vs n_probes")
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------
# TEST 4: Cholesky SE with A_new (linear combinations) vs coefficient SE
# -----------------------------------------------------------------------
cat("\n--- Test 4: posterior_se with A_new ---\n")

n <- 2000; p <- 500
prob  <- make_problem(n, p)
fit_c <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")

n_new_sizes <- c(10, 100, 500, 1000, 5000)
results_lc <- data.frame()

for (n_new in n_new_sizes) {
  A_new <- matrix(rnorm(n_new * p), n_new, p)
  t <- system.time({
    se_lc <- posterior_se(fit_c, A_new = A_new)
  })["elapsed"]
  cat(sprintf("  n_new = %5d: %.2fs\n", n_new, t))
  results_lc <- rbind(results_lc, data.frame(n_new = n_new, time = t))
}

plot(results_lc$n_new, results_lc$time,
     type = "b", pch = 19, col = "#762a83",
     xlab = "n_new (rows of A_new)", ylab = "Time (s)",
     main = "posterior_se with A_new: time vs n_new\n(Cholesky path, p = 500)")

cat("\n=== Speed tests complete ===\n")
dev.off()
cat("Plots saved to vignettes/speed_tests.pdf\n")
