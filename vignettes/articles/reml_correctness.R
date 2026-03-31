# REML correctness test
# Use models where the marginal likelihood has a known analytic form
# Compare tune_reml output to exact analytic optimum
# Run with: pkgload::load_all("."); source("vignettes/reml_correctness.R")

pkgload::load_all(".")
library(Matrix)

pdf("vignettes/reml_correctness.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

cat("=== REML correctness tests ===\n\n")

# -----------------------------------------------------------------------
# TEST 1: A = I, Q = I (ridge prior)
# y ~ N(0, sigma2e * (1 + phi) * I)
# Marginal ll: -n/2 * log(sigma2e*(1+phi)) - y'y / (2*sigma2e*(1+phi))
# Analytic optimum: sigma2e*(1+phi) = y'y/n  (MLE of variance)
# phi is not separately identified from sigma2e -- only phi*sigma2e is
# So we fix sigma2e=1 and estimate phi
# -----------------------------------------------------------------------
cat("--- Test 1: A=I, Q=I, analytic marginal ll ---\n")

set.seed(42)
n       <- 200
p       <- 200   # A = I so n = p
phi_true    <- 4
sigma2e_true <- 1

# y ~ N(0, sigma2e*(1+phi)*I)
y1 <- rnorm(n, sd = sqrt(sigma2e_true * (1 + phi_true)))
A1 <- diag(n)
Q1 <- Matrix::Diagonal(n, 1)

# Analytic marginal ll as function of phi (with sigma2e profiled)
analytic_ll <- function(phi, y) {
  n    <- length(y)
  s2   <- sum(y^2) / (n * (1 + phi))   # profiled sigma2e
  if (s2 <= 0) return(-Inf)
  -n/2 * log(s2) - n/2 * log(1 + phi) - n/2
}

phi_grid    <- seq(0.5, 15, by = 0.5)
ll_analytic <- sapply(phi_grid, analytic_ll, y = y1)
phi_analytic_opt <- phi_grid[which.max(ll_analytic)]
# exact analytic optimum: phi = var(y)/sigma2e - 1
phi_exact <- var(y1) - 1
cat(sprintf("  Analytic optimum phi: %.3f\n", phi_exact))

# Q_fun for tune_reml
Q_fun1 <- function(theta) list(Q = Q1, Q_matrix = Q1)

tuned1 <- tune_reml(
  y1, A1, Q_fun1,
  theta_init    = numeric(0),
  logdet_method = "cholesky",
  verbose       = FALSE
)
cat(sprintf("  tune_reml phi:        %.3f\n", tuned1$phi))
cat(sprintf("  Difference:           %.4f\n\n", abs(tuned1$phi - phi_exact)))

# plot
plot(phi_grid, ll_analytic - max(ll_analytic),
     type="b", pch=19, col="#2166ac",
     xlab="phi", ylab="ll (normalized)",
     main="Test 1: A=I, Q=I\nAnalytic vs tune_reml")
abline(v=phi_exact,    col="red",    lty=2, lwd=2)
abline(v=tuned1$phi,   col="#4dac26", lty=2, lwd=2)
legend("bottomleft",
       c(sprintf("Analytic opt=%.2f", phi_exact),
         sprintf("tune_reml=%.2f",   tuned1$phi)),
       col=c("red","#4dac26"), lty=2, bty="n")

# -----------------------------------------------------------------------
# TEST 2: A = I, Q = tau*I, estimate tau (shape of prior)
# y ~ N(0, sigma2e*(1 + phi/tau)*I)
# Now tau is in theta -- should be identifiable from phi/tau ratio
# -----------------------------------------------------------------------
cat("--- Test 2: A=I, Q=tau*I, estimate tau ---\n")

tau_true <- 2
y2 <- rnorm(n, sd = sqrt(sigma2e_true * (1 + phi_true/tau_true)))

analytic_ll_tau <- function(tau, phi_range, y) {
  # profile phi at each tau
  best_ll <- -Inf
  for (phi in phi_range) {
    s2 <- sum(y^2) / (n * (1 + phi/tau))
    if (s2 <= 0) next
    ll <- -n/2*log(s2) - n/2*log(1+phi/tau) - n/2 +
      n/2*log(tau) - n/2*log(phi)   # logdetQ = n*log(tau), -p/2*log(phi)
    if (ll > best_ll) best_ll <- ll
  }
  best_ll
}

tau_grid    <- seq(0.5, 8, by=0.25)
phi_range   <- seq(0.1, 50, by=0.1)
ll_tau      <- sapply(tau_grid, analytic_ll_tau, phi_range=phi_range, y=y2)
tau_analytic_opt <- tau_grid[which.max(ll_tau)]
cat(sprintf("  True tau: %.1f,  Analytic optimum tau: %.3f\n", tau_true, tau_analytic_opt))

Q_fun2 <- function(theta) {
  tau <- theta[["tau"]]
  if (tau <= 0) return(list(Q=NULL))
  Q <- Matrix::Diagonal(n, tau)
  list(Q=Q, Q_matrix=Q,
       log_det_Q = n * log(tau))
}

tuned2 <- tune_reml(
  y2, A1, Q_fun2,
  theta_init    = c(tau=1),
  lower         = 0.01,
  upper         = 20,
  logdet_method = "cholesky",
  verbose       = FALSE
)
cat(sprintf("  tune_reml tau: %.3f\n", tuned2$theta[["tau"]]))
cat(sprintf("  Difference:    %.4f\n\n", abs(tuned2$theta[["tau"]] - tau_true)))

plot(tau_grid, ll_tau - max(ll_tau),
     type="b", pch=19, col="#2166ac",
     xlab="tau", ylab="ll (normalized)",
     main="Test 2: A=I, Q=tau*I\nAnalytic vs tune_reml")
abline(v=tau_true,                   col="red",     lty=2, lwd=2)
abline(v=tau_analytic_opt,           col="gray",    lty=2, lwd=2)
abline(v=tuned2$theta[["tau"]],      col="#4dac26", lty=2, lwd=2)
legend("bottomleft",
       c(sprintf("True=%.1f", tau_true),
         sprintf("Analytic=%.2f", tau_analytic_opt),
         sprintf("tune_reml=%.2f", tuned2$theta[["tau"]])),
       col=c("red","gray","#4dac26"), lty=2, bty="n")

# -----------------------------------------------------------------------
# TEST 3: general A, Q=I -- estimate phi only
# Marginal ll: y ~ N(0, sigma2e*(I + phi*A*A'))
# This is the standard ridge regression marginal likelihood
# -----------------------------------------------------------------------
cat("--- Test 3: general A, Q=I, estimate phi ---\n")

set.seed(42)
n3   <- 100
p3   <- 50
phi3 <- 3
A3   <- matrix(rnorm(n3*p3), n3, p3) / sqrt(p3)
Q3   <- Matrix::Diagonal(p3, 1)

x3   <- rnorm(p3, sd=sqrt(phi3))
y3   <- as.numeric(A3 %*% x3 + rnorm(n3))

# Exact marginal ll via direct computation
exact_ll_phi <- function(phi, y, A) {
  n    <- length(y)
  Sy   <- diag(n) + phi * A %*% t(A)
  ld   <- as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
  # profile sigma2e
  s2   <- as.numeric(t(y) %*% solve(Sy, y)) / n
  -n/2*log(s2) - 1/2*ld - n/2*log(s2) + 1/2*ld   # simplifies
  # actually: -n/2*log(2*pi*s2) - 1/2*ld - n/(2*s2) * y'Sy^{-1}y / s2...
  # let's just compute it directly
  -n/2*log(s2) - 1/2*ld - n/2
}

# recompute cleanly
exact_ll_phi <- function(phi, y, A) {
  n  <- length(y)
  Sy <- diag(n) + phi * tcrossprod(A)
  ld <- as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
  s2 <- as.numeric(t(y) %*% solve(Sy, y)) / n
  if (s2 <= 0) return(-Inf)
  -n/2*log(s2) - 1/2*ld - n/2
}

phi_grid3   <- seq(0.2, 10, by=0.2)
ll_exact3   <- sapply(phi_grid3, exact_ll_phi, y=y3, A=A3)
phi_exact3  <- phi_grid3[which.max(ll_exact3)]
cat(sprintf("  True phi: %.1f,  Exact optimum phi: %.3f\n", phi3, phi_exact3))

Q_fun3 <- function(theta) list(Q=Q3, Q_matrix=Q3)
tuned3 <- tune_reml(
  y3, A3, Q_fun3,
  theta_init    = numeric(0),
  logdet_method = "cholesky",
  verbose       = FALSE
)
cat(sprintf("  tune_reml phi: %.3f\n", tuned3$phi))
cat(sprintf("  Difference:    %.4f\n\n", abs(tuned3$phi - phi_exact3)))

plot(phi_grid3, ll_exact3 - max(ll_exact3),
     type="b", pch=19, col="#2166ac",
     xlab="phi", ylab="ll (normalized)",
     main="Test 3: general A, Q=I\nExact marginal ll vs tune_reml")
abline(v=phi3,         col="red",     lty=2, lwd=2)
abline(v=phi_exact3,   col="gray",    lty=2, lwd=2)
abline(v=tuned3$phi,   col="#4dac26", lty=2, lwd=2)
legend("bottomleft",
       c(sprintf("True=%.1f", phi3),
         sprintf("Exact=%.2f", phi_exact3),
         sprintf("tune_reml=%.2f", tuned3$phi)),
       col=c("red","gray","#4dac26"), lty=2, bty="n")

# -----------------------------------------------------------------------
# TEST 4: scan over reml ll internally vs exact, point by point
# Most direct correctness check: for each phi, does tune_reml's
# internal ll match the exact marginal ll?
# -----------------------------------------------------------------------
cat("--- Test 4: point-by-point ll comparison ---\n")

# compute tune_reml's internal ll at each phi_grid3 value
Rinv        <- resolve_Rinv(NULL, n3)
AtRinvA     <- Matrix::crossprod(A3)
apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
AtRinvy     <- as.numeric(Matrix::crossprod(A3, y3))
yRinvy      <- as.numeric(crossprod(y3, y3))
probes      <- matrix(sample(c(-1L,1L), p3*50, replace=TRUE), p3, 50)
logdetQ     <- 0   # Q=I so log|Q|=0

apply_Q3 <- as_apply(Q3)
prior3   <- list(Q=Q3, Q_matrix=Q3, AtRinvA_matrix=AtRinvA)

ll_internal <- sapply(phi_grid3, function(phi) {
  res <- fastblm:::.eval_reml_ll(
    phi, AtRinvy, yRinvy, apply_AtRinvA, apply_Q3, prior3,
    p3, n3, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p3, rep(0,p3)
  )
  res$ll
})

# plot both -- normalize to same reference (exact ll at phi=1)
ref_exact   <- ll_exact[which(phi_grid3 == 1)]
ref_internal <- ll_internal[which(phi_grid3 == 1)]

plot(phi_grid3, ll_exact3 - ref_exact,
     type="b", pch=19, col="#2166ac",
     xlab="phi", ylab="ll (normalized to phi=1)",
     main="Test 4: exact vs internal ll\npoint by point (same normalization)")
lines(phi_grid3, ll_internal - ref_internal,
      type="b", pch=19, col="#e66101")
abline(v=phi_exact3, col="red", lty=2)
legend("bottomleft",
       c("Exact marginal ll", "tune_reml internal ll"),
       col=c("#2166ac","#e66101"), pch=19, lty=1, bty="n")

cat(sprintf("  Exact opt phi:    %.2f\n", phi_exact3))
cat(sprintf("  Internal opt phi: %.2f\n", phi_grid3[which.max(ll_internal)]))
cat(sprintf("  Constant offset (should be -n/2*log(2pi) - n/2 = %.4f): %.4f\n",
            -n3/2*log(2*pi) - n3/2,
            mean(ll_exact3 - ll_internal)))

# -----------------------------------------------------------------------
# TEST 5: sigma2e recovery
# At the true phi, does profiled sigma2e match truth?
# -----------------------------------------------------------------------
cat("--- Test 5: sigma2e recovery ---\n")

n5 <- 500; p5 <- 100
set.seed(42)
phi5     <- 4
sigma2e5 <- 2

A5   <- matrix(rnorm(n5*p5), n5, p5) / sqrt(p5)
Q5   <- Matrix::Diagonal(p5, 1)
x5   <- rnorm(p5, sd=sqrt(phi5*sigma2e5))
y5   <- as.numeric(A5 %*% x5 + rnorm(n5, sd=sqrt(sigma2e5)))

Q_fun5 <- function(theta) list(Q=Q5, Q_matrix=Q5)
tuned5 <- tune_reml(y5, A5, Q_fun5,
                    theta_init=numeric(0),
                    logdet_method="cholesky",
                    verbose=FALSE)

cat(sprintf("  True phi=%.1f  estimated=%.3f\n", phi5, tuned5$phi))
cat(sprintf("  True sigma2e=%.1f  estimated=%.3f\n", sigma2e5, tuned5$sigma2e))
cat(sprintf("  True sigma2b=%.1f  estimated=%.3f\n",
            phi5*sigma2e5, tuned5$sigma2b))

# verify sigma2e from fit matches tune
fit5 <- fit_fastblm(y5, A5, Q5, phi=tuned5$phi, solver="cholesky")
cat(sprintf("  sigma2e from fit: %.3f  from tune: %.3f  diff: %.4f\n\n",
            fit5$sigma2e, tuned5$sigma2e, abs(fit5$sigma2e - tuned5$sigma2e)))

# plot phi surface with sigma2e shown
phi_grid5 <- seq(1, 12, by=0.25)
AtRinvA5  <- Matrix::crossprod(A5)
apply_AA5 <- function(v) as.numeric(AtRinvA5 %*% v)
AtRinvy5  <- as.numeric(Matrix::crossprod(A5, y5))
yRinvy5   <- as.numeric(crossprod(y5))
prior5    <- list(Q=Q5, Q_matrix=Q5, AtRinvA_matrix=AtRinvA5)
apply_Q5  <- as_apply(Q5)
probes5   <- matrix(sample(c(-1L,1L), p5*50, replace=TRUE), p5, 50)

ll5 <- sapply(phi_grid5, function(phi) {
  res <- fastblm:::.eval_reml_ll(
    phi, AtRinvy5, yRinvy5, apply_AA5, apply_Q5, prior5,
    p5, n5, 0, probes5, 50L, "cholesky", 1e-6, 4L*p5, rep(0,p5)
  )
  res$ll
})

plot(phi_grid5, ll5 - max(ll5),
     type="b", pch=19, col="#2166ac",
     xlab="phi", ylab="ll (normalized)",
     main=sprintf("Test 5: sigma2e recovery\ntrue phi=%.1f, sigma2e=%.1f", phi5, sigma2e5))
abline(v=phi5,      col="red",     lty=2, lwd=2)
abline(v=tuned5$phi, col="#4dac26", lty=2, lwd=2)
legend("bottomleft",
       c(sprintf("True phi=%.1f", phi5),
         sprintf("Estimated phi=%.2f", tuned5$phi)),
       col=c("red","#4dac26"), lty=2, bty="n")

par(mfrow=c(1,1))
dev.off()
cat("=== REML correctness tests complete ===\n")
cat("Plots saved to vignettes/reml_correctness.pdf\n")
