# Targeted debug: check quadratic form and logdet separately
pkgload::load_all(".")
library(Matrix)

set.seed(42)
n <- 50; p <- 20
A <- matrix(rnorm(n*p), n, p) / sqrt(p)
Q <- Matrix::Diagonal(p, 1)
y <- rnorm(n, sd=2)

cat("=== Term by term comparison ===\n\n")

for (phi in c(0.5, 1, 2, 3, 5, 8)) {
  # Exact quantities in observation space
  Sy       <- diag(n) + phi * tcrossprod(A)
  Sy_inv   <- solve(Sy)
  qform_exact <- as.numeric(t(y) %*% Sy_inv %*% y)
  logdet_Sy   <- as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
  s2_exact    <- qform_exact / n

  # Our quantities in coefficient space
  K        <- t(A) %*% A + (1/phi) * diag(p)
  Aty      <- as.numeric(t(A) %*% y)
  x_hat    <- solve(K, Aty)
  qform_ours  <- sum(y^2) - as.numeric(t(Aty) %*% x_hat)
  logdet_K    <- as.numeric(determinant(K, logarithm=TRUE)$modulus)
  s2_ours     <- qform_ours / n

  # What logdet_Sy should equal via det lemma
  logdet_Sy_via_K <- p*log(phi) + logdet_K  # since Q=I, logdetQ=0

  cat(sprintf("phi=%.1f:\n", phi))
  cat(sprintf("  qform: exact=%.4f  ours=%.4f  match=%s\n",
              qform_exact, qform_ours, isTRUE(all.equal(qform_exact, qform_ours))))
  cat(sprintf("  logdet_Sy: exact=%.4f  via_K=%.4f  match=%s\n",
              logdet_Sy, logdet_Sy_via_K,
              isTRUE(all.equal(logdet_Sy, logdet_Sy_via_K))))
  cat(sprintf("  s2: exact=%.4f  ours=%.4f\n", s2_exact, s2_ours))

  # Exact profiled ll (without 2pi constant)
  ll_exact_nopi <- -n/2*log(s2_exact) - 1/2*logdet_Sy - n/2
  # Our ll
  ll_ours <- -n/2*log(s2_ours) - 1/2*logdet_K - p/2*log(phi) + 0  # logdetQ=0
  # Should differ by exactly -n/2*log(2pi)
  diff    <- ll_exact_nopi - ll_ours
  cat(sprintf("  ll_exact(no2pi)=%.4f  ll_ours=%.4f  diff=%.4f  expected=%.4f\n\n",
              ll_exact_nopi, ll_ours, diff, -n/2))
}
