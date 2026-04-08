# Hutchinson diagonal estimator test
# Tests speed and accuracy of Hutchinson SE vs Cholesky SE
# Run with: pkgload::load_all("."); source("vignettes/hutchinson_test.R")

pkgload::load_all(".")
library(Matrix)

cat("=== Hutchinson SE test ===\n\n")

set.seed(42)
n <- 1000
p <- 200

# simple dense A, diagonal Q
A <- matrix(rnorm(n * p), n, p)
Q <- Matrix::Diagonal(p, 1)
phi <- 5

# fit both ways
y     <- as.numeric(A %*% rnorm(p) + rnorm(n))
fit_c <- fit_fastblm(y, A, Q, phi = phi, solver = "cholesky")
fit_p <- fit_fastblm(y, A, Q, phi = phi, solver = "pcg")

# exact SE from Cholesky
se_exact <- posterior_se(fit_c)
cat("Exact SE range:", range(se_exact), "\n\n")

# test Hutchinson for increasing n_probes
cat("n_probes | time(s) | cor  | mean_rel_err\n")
cat("---------+---------+------+-------------\n")
for (np in c(10, 25, 50, 100, 200, 500, 1000)) {
  t <- system.time(se_h <- posterior_se(fit_p, n_probes = np))["elapsed"]
  cat(sprintf("%8d | %7.2f | %.3f | %.4f\n",
              np, t, cor(se_h, se_exact), mean(abs(se_h - se_exact)/se_exact)))
}
