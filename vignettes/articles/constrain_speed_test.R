# constrain() speed test
# Key question: is constrain() fast when q (number of constraints) is small?
# Tests across problem sizes and numbers of constraints
# Run with: pkgload::load_all("."); source("vignettes/constrain_speed_test.R")

pkgload::load_all(".")
library(Matrix)

pdf("vignettes/constrain_speed_test.pdf", width = 12, height = 8)

cat("=== constrain() speed test ===\n\n")

# -----------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------
make_problem <- function(n, p, phi = 5, rho = 0.5, seed = 42) {
  set.seed(seed)
  A <- matrix(rnorm(n * p), n, p)
  Q <- Matrix::bandSparse(p, k = c(-1, 0, 1),
                          diagonals = list(rep(-rho, p-1),
                                           c(1, rep(1 + rho^2, p-2), 1),
                                           rep(-rho, p-1)))
  Q <- Matrix::forceSymmetric(Q)
  y <- rnorm(n)
  list(y = y, A = A, Q = Q, phi = phi)
}

make_C <- function(q, p) {
  # q random linear constraints
  C <- matrix(rnorm(q * p), q, p)
  C
}

# -----------------------------------------------------------------------
# TEST 1: constrain speed vs problem size p, fixed q = 1
# -----------------------------------------------------------------------
cat("--- Test 1: constrain speed vs p, q = 1 ---\n")

p_sizes <- c(100, 500, 1000, 2000, 5000)
results1 <- data.frame()

for (p in p_sizes) {
  n    <- max(50, p %/% 5)
  prob <- make_problem(n, p)
  fit  <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")
  C    <- make_C(1, p)

  t_con <- system.time(fit_c <- constrain(fit, C))["elapsed"]
  t_se  <- system.time(posterior_se(fit_c))["elapsed"]

  cat(sprintf("  p = %5d, n = %5d: constrain = %.3fs  posterior_se = %.3fs\n",
              p, n, t_con, t_se))

  results1 <- rbind(results1, data.frame(p=p, n=n, q=1, t_con=t_con, t_se=t_se))
}

# -----------------------------------------------------------------------
# TEST 2: constrain speed vs q (number of constraints), fixed p
# -----------------------------------------------------------------------
cat("\n--- Test 2: constrain speed vs q, fixed p = 1000 ---\n")

p_fixed <- 1000
n_fixed <- 200
prob    <- make_problem(n_fixed, p_fixed)
fit     <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")

q_sizes <- c(1, 2, 5, 10, 20, 50)
results2 <- data.frame()

for (q in q_sizes) {
  C     <- make_C(q, p_fixed)
  t_con <- system.time(fit_c <- constrain(fit, C))["elapsed"]
  t_se  <- system.time(posterior_se(fit_c))["elapsed"]

  cat(sprintf("  q = %3d: constrain = %.3fs  posterior_se = %.3fs\n", q, t_con, t_se))
  results2 <- rbind(results2, data.frame(q=q, t_con=t_con, t_se=t_se))
}

# -----------------------------------------------------------------------
# TEST 3: constrain + posterior_se vs p, realistic q = 1 and q = 2
# comparing to unconstrained posterior_se
# -----------------------------------------------------------------------
cat("\n--- Test 3: constrained vs unconstrained posterior_se ---\n")

results3 <- data.frame()

for (p in c(100, 500, 1000, 2000)) {
  n    <- max(50, p %/% 5)
  prob <- make_problem(n, p)
  fit  <- fit_fastblm(prob$y, prob$A, prob$Q, phi = prob$phi, solver = "cholesky")

  t_unc  <- system.time(posterior_se(fit))["elapsed"]

  C1     <- make_C(1, p)
  fit_c1 <- constrain(fit, C1)
  t_con1 <- system.time(posterior_se(fit_c1))["elapsed"]

  C2     <- make_C(2, p)
  fit_c2 <- constrain(fit, C2)
  t_con2 <- system.time(posterior_se(fit_c2))["elapsed"]

  cat(sprintf("  p = %5d: unconstrained = %.3fs  q=1 = %.3fs  q=2 = %.3fs\n",
              p, t_unc, t_con1, t_con2))

  results3 <- rbind(results3, data.frame(
    p=p, t_unc=t_unc, t_con1=t_con1, t_con2=t_con2
  ))
}

# -----------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------
par(mfrow = c(1, 3))

# Plot 1: constrain time vs p
plot(results1$p, results1$t_con,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = "constrain() time vs p\n(q = 1)",
     ylim = range(c(results1$t_con, results1$t_se)))
lines(results1$p, results1$t_se, type = "b", pch = 19, col = "#e66101")
legend("topleft", legend = c("constrain()", "posterior_se()"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

# Plot 2: time vs q
plot(results2$q, results2$t_con,
     type = "b", pch = 19, col = "#2166ac",
     xlab = "q (number of constraints)", ylab = "Time (s)",
     main = sprintf("constrain() time vs q\n(p = %d)", p_fixed),
     ylim = range(c(results2$t_con, results2$t_se)))
lines(results2$q, results2$t_se, type = "b", pch = 19, col = "#e66101")
legend("topleft", legend = c("constrain()", "posterior_se()"),
       col = c("#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

# Plot 3: constrained vs unconstrained SE time
plot(results3$p, results3$t_unc,
     type = "b", pch = 19, col = "#4dac26",
     xlab = "p (coefficients)", ylab = "Time (s)",
     main = "posterior_se: constrained vs unconstrained",
     ylim = range(c(results3$t_unc, results3$t_con1, results3$t_con2)))
lines(results3$p, results3$t_con1, type = "b", pch = 19, col = "#2166ac")
lines(results3$p, results3$t_con2, type = "b", pch = 19, col = "#e66101")
legend("topleft", legend = c("unconstrained", "q=1", "q=2"),
       col = c("#4dac26", "#2166ac", "#e66101"), pch = 19, lty = 1, bty = "n")

par(mfrow = c(1, 1))

cat("\n=== constrain speed test complete ===\n")
dev.off()
cat("Plots saved to vignettes/constrain_speed_test.pdf\n")
