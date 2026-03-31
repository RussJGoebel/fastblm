# SAR identifiability: does n >= p help?
# Key question: is rho identifiable when n >= p with random A?
# Prediction: no, because random A averages out spatial structure
# Run with: pkgload::load_all("."); source("vignettes/sar_identifiability_np.R")

pkgload::load_all(".")
library(Matrix)

make_SAR_Q <- function(p, rho) {
  i_idx    <- c(2:p, 1:(p-1))
  j_idx    <- c(1:(p-1), 2:p)
  row_sums <- c(1, rep(2, p-2), 1)
  x_vals   <- 1 / row_sums[i_idx]
  W <- Matrix::sparseMatrix(i=i_idx, j=j_idx, x=x_vals, dims=c(p,p))
  S <- Matrix::Diagonal(p) - rho * W
  Matrix::forceSymmetric(Matrix::crossprod(S))
}

p        <- 50
rho_true <- 0.7
phi_true <- 5
sigma2e_true <- 0.5
rho_grid <- seq(0.05, 0.95, by = 0.05)

set.seed(42)
Q_true <- make_SAR_Q(p, rho_true)
x_true <- as.numeric(Matrix::solve(Matrix::Cholesky(Q_true), rnorm(p))) *
  sqrt(phi_true * sigma2e_true)

pdf("vignettes/sar_identifiability_np.pdf", width = 14, height = 10)
par(mfrow = c(2, 4))

# -----------------------------------------------------------------------
# TEST 1: random A, varying n/p ratio
# -----------------------------------------------------------------------
cat("--- Test 1: random A, varying n/p ratio ---\n")

n_sizes <- c(25, 50, 100, 500)  # n/p = 0.5, 1, 2, 10

for (n in n_sizes) {
  set.seed(42)
  A <- matrix(rnorm(n * p), n, p) / sqrt(p)
  y <- as.numeric(A %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

  AtRinvA       <- Matrix::crossprod(A)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
  AtRinvy       <- as.numeric(Matrix::crossprod(A, y))
  yRinvy        <- as.numeric(crossprod(y, y))
  probes        <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)

  ll_vals <- sapply(rho_grid, function(rho) {
    Q       <- make_SAR_Q(p, rho)
    apply_Q <- as_apply(Q)
    logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
    prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA)

    phi_i <- fastblm:::.profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, log(0.01), log(1000),
      "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res <- fastblm:::.eval_reml_ll(
      phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res$ll
  })

  peak_rho <- rho_grid[which.max(ll_vals)]
  cat(sprintf("  n/p = %.1f: peak rho = %.2f\n", n/p, peak_rho))

  plot(rho_grid, ll_vals, type="b", pch=19, col="#2166ac",
       xlab="rho", ylab="REML ll",
       main=sprintf("Random A, n/p = %.1f\npeak rho = %.2f", n/p, peak_rho))
  abline(v=rho_true, col="red", lty=2)
  abline(v=peak_rho, col="#4dac26", lty=2)
  legend("bottomright",
         c(sprintf("True=%.1f", rho_true), sprintf("Peak=%.2f", peak_rho)),
         col=c("red","#4dac26"), lty=2, bty="n", cex=0.8)
}

# -----------------------------------------------------------------------
# TEST 2: identity A (A = I, so y = x + epsilon directly)
# This is maximum alignment -- should be most identifiable
# -----------------------------------------------------------------------
cat("\n--- Test 2: A = I (direct observation) ---\n")
n <- p
A_ident <- diag(p)
y_ident <- as.numeric(A_ident %*% x_true + rnorm(p, sd=sqrt(sigma2e_true)))

AtRinvA_ident   <- Matrix::crossprod(A_ident)
apply_AtRinvA_i <- function(v) as.numeric(AtRinvA_ident %*% v)
AtRinvy_ident   <- as.numeric(Matrix::crossprod(A_ident, y_ident))
yRinvy_ident    <- as.numeric(crossprod(y_ident, y_ident))

ll_ident <- sapply(rho_grid, function(rho) {
  Q       <- make_SAR_Q(p, rho)
  apply_Q <- as_apply(Q)
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
  prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA_ident)

  phi_i <- fastblm:::.profile_phi(
    AtRinvy_ident, yRinvy_ident, apply_AtRinvA_i, apply_Q, prior,
    p, p, logdetQ, probes, 50L, log(0.01), log(1000),
    "cholesky", 1e-6, 4L*p, rep(0,p)
  )
  res <- fastblm:::.eval_reml_ll(
    phi_i, AtRinvy_ident, yRinvy_ident, apply_AtRinvA_i, apply_Q, prior,
    p, p, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
  )
  res$ll
})

peak_rho_ident <- rho_grid[which.max(ll_ident)]
cat(sprintf("  A=I: peak rho = %.2f (true = %.2f)\n", peak_rho_ident, rho_true))

plot(rho_grid, ll_ident, type="b", pch=19, col="#762a83",
     xlab="rho", ylab="REML ll",
     main=sprintf("A = I (direct observation)\npeak rho = %.2f", peak_rho_ident))
abline(v=rho_true, col="red", lty=2)
abline(v=peak_rho_ident, col="#4dac26", lty=2)
legend("bottomright",
       c(sprintf("True=%.1f", rho_true), sprintf("Peak=%.2f", peak_rho_ident)),
       col=c("red","#4dac26"), lty=2, bty="n", cex=0.8)

# -----------------------------------------------------------------------
# TEST 3: smooth spatial A (each obs = weighted average of neighbors)
# More structured than random, less perfect than identity
# -----------------------------------------------------------------------
cat("\n--- Test 3: smooth spatial A ---\n")
n <- 200
set.seed(42)
# each obs is a gaussian-weighted average of nearby basis functions
A_smooth <- matrix(0, n, p)
obs_locs  <- seq(0, 1, length.out=n)
bas_locs  <- seq(0, 1, length.out=p)
bw        <- 3/p  # bandwidth
for (i in seq_len(n)) {
  w <- exp(-0.5 * ((obs_locs[i] - bas_locs)/bw)^2)
  A_smooth[i,] <- w / sum(w)
}

y_smooth <- as.numeric(A_smooth %*% x_true + rnorm(n, sd=sqrt(sigma2e_true)))

AtRinvA_sm    <- Matrix::crossprod(A_smooth)
apply_AtRinvA_sm <- function(v) as.numeric(AtRinvA_sm %*% v)
AtRinvy_sm    <- as.numeric(Matrix::crossprod(A_smooth, y_smooth))
yRinvy_sm     <- as.numeric(crossprod(y_smooth, y_smooth))

ll_smooth <- sapply(rho_grid, function(rho) {
  Q       <- make_SAR_Q(p, rho)
  apply_Q <- as_apply(Q)
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
  prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA_sm)

  phi_i <- fastblm:::.profile_phi(
    AtRinvy_sm, yRinvy_sm, apply_AtRinvA_sm, apply_Q, prior,
    p, n, logdetQ, probes, 50L, log(0.01), log(1000),
    "cholesky", 1e-6, 4L*p, rep(0,p)
  )
  res <- fastblm:::.eval_reml_ll(
    phi_i, AtRinvy_sm, yRinvy_sm, apply_AtRinvA_sm, apply_Q, prior,
    p, n, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
  )
  res$ll
})

peak_rho_smooth <- rho_grid[which.max(ll_smooth)]
cat(sprintf("  Smooth A: peak rho = %.2f (true = %.2f)\n", peak_rho_smooth, rho_true))

plot(rho_grid, ll_smooth, type="b", pch=19, col="#e66101",
     xlab="rho", ylab="REML ll",
     main=sprintf("Smooth spatial A, n=%d\npeak rho = %.2f", n, peak_rho_smooth))
abline(v=rho_true, col="red", lty=2)
abline(v=peak_rho_smooth, col="#4dac26", lty=2)
legend("bottomright",
       c(sprintf("True=%.1f", rho_true), sprintf("Peak=%.2f", peak_rho_smooth)),
       col=c("red","#4dac26"), lty=2, bty="n", cex=0.8)

par(mfrow=c(1,1))
dev.off()
cat("Plots saved to vignettes/sar_identifiability_np.pdf\n")
