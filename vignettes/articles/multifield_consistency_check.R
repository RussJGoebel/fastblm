# Consistency check: multiple independent field realizations
# Key insight: with one field realization, rho is identified from a single
# sample path. With multiple realizations, we have genuine replication.
# Run with: pkgload::load_all("."); source("vignettes/multi_field_consistency.R")

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

make_smooth_A <- function(n, p, bw = 3/p) {
  obs_locs <- seq(0, 1, length.out = n)
  bas_locs <- seq(0, 1, length.out = p)
  A <- matrix(0, n, p)
  for (i in seq_len(n)) {
    w <- exp(-0.5 * ((obs_locs[i] - bas_locs)/bw)^2)
    A[i,] <- w / sum(w)
  }
  A
}

p            <- 50
rho_true     <- 0.7
phi_true     <- 5
sigma2e_true <- 0.5
rho_grid     <- seq(0.05, 0.90, by = 0.05)

Q_true <- make_SAR_Q(p, rho_true)
CQ     <- Matrix::Cholesky(Q_true)
probes <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)

# observations per field realization
n_per_field <- 50
A <- make_smooth_A(n_per_field, p)

cat("=== Multi-field consistency check ===\n\n")
cat(sprintf("p = %d, n_per_field = %d, true rho = %.1f\n\n", p, n_per_field, rho_true))

pdf("vignettes/multi_field_consistency.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# -----------------------------------------------------------------------
# TEST 1: ll surface vs number of field realizations
# Stack observations from multiple fields: y = [A x_1; A x_2; ...]
# -----------------------------------------------------------------------
cat("--- Test 1: ll surface vs number of field realizations ---\n")

n_fields_list <- c(1, 5, 20, 100, 500)
cols <- c("#2166ac","#4dac26","#e66101","#762a83","#d01c8b")

ll_all <- matrix(NA, length(n_fields_list), length(rho_grid))

for (fi in seq_along(n_fields_list)) {
  n_fields <- n_fields_list[fi]
  set.seed(42)

  # stack n_fields independent realizations
  y_list <- lapply(seq_len(n_fields), function(k) {
    x_k <- as.numeric(Matrix::solve(CQ, rnorm(p))) * sqrt(phi_true * sigma2e_true)
    as.numeric(A %*% x_k + rnorm(n_per_field, sd=sqrt(sigma2e_true)))
  })

  # block diagonal A: same A repeated n_fields times
  A_stack <- do.call(rbind, replicate(n_fields, A, simplify=FALSE))
  y_stack <- unlist(y_list)
  n_total <- nrow(A_stack)

  AtRinvA       <- Matrix::crossprod(A_stack)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
  AtRinvy       <- as.numeric(Matrix::crossprod(A_stack, y_stack))
  yRinvy        <- as.numeric(crossprod(y_stack, y_stack))

  ll_vals <- sapply(rho_grid, function(rho) {
    Q       <- make_SAR_Q(p, rho)
    apply_Q <- as_apply(Q)
    logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
    prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA)

    phi_i <- fastblm:::.profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n_total, logdetQ, probes, 50L, log(0.01), log(1000),
      "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res <- fastblm:::.eval_reml_ll(
      phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n_total, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
    )
    res$ll
  })

  ll_all[fi,] <- ll_vals - max(ll_vals)
  peak_rho <- rho_grid[which.max(ll_vals)]
  cat(sprintf("  n_fields = %3d (n_total = %5d): peak rho = %.2f\n",
              n_fields, n_total, peak_rho))
}

plot(rho_grid, ll_all[1,], type="b", pch=19, col=cols[1],
     xlab="rho", ylab="REML ll (normalized)",
     main="ll surface vs number of field realizations\n(normalized to peak = 0)",
     ylim=range(ll_all, na.rm=TRUE))
for (fi in 2:length(n_fields_list))
  lines(rho_grid, ll_all[fi,], type="b", pch=19, col=cols[fi])
abline(v=rho_true, col="red", lty=2)
legend("bottomleft", sprintf("%d fields", n_fields_list),
       col=cols, pch=19, lty=1, bty="n", cex=0.8)

# -----------------------------------------------------------------------
# TEST 2: replicated estimation with multiple fields
# -----------------------------------------------------------------------
cat("\n--- Test 2: replicated estimation ---\n")

n_fields_sizes <- c(1, 5, 20, 100, 500)
n_reps         <- 20L
results        <- data.frame()

for (n_fields in n_fields_sizes) {
  rho_ests <- numeric(n_reps)
  n_total  <- n_fields * n_per_field

  for (rep in seq_len(n_reps)) {
    set.seed(rep * 100)

    y_stack <- unlist(lapply(seq_len(n_fields), function(k) {
      x_k <- as.numeric(Matrix::solve(CQ, rnorm(p))) * sqrt(phi_true * sigma2e_true)
      as.numeric(A %*% x_k + rnorm(n_per_field, sd=sqrt(sigma2e_true)))
    }))

    A_stack       <- do.call(rbind, replicate(n_fields, A, simplify=FALSE))
    AtRinvA       <- Matrix::crossprod(A_stack)
    apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)
    AtRinvy       <- as.numeric(Matrix::crossprod(A_stack, y_stack))
    yRinvy        <- as.numeric(crossprod(y_stack, y_stack))

    ll_vals <- sapply(rho_grid, function(rho) {
      Q       <- make_SAR_Q(p, rho)
      apply_Q <- as_apply(Q)
      logdetQ <- as.numeric(Matrix::determinant(Q, logarithm=TRUE)$modulus)
      prior   <- list(Q=Q, Q_matrix=Q, AtRinvA_matrix=AtRinvA)

      phi_i <- fastblm:::.profile_phi(
        AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n_total, logdetQ, probes, 50L, log(0.01), log(1000),
        "cholesky", 1e-6, 4L*p, rep(0,p)
      )
      res <- fastblm:::.eval_reml_ll(
        phi_i, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n_total, logdetQ, probes, 50L, "cholesky", 1e-6, 4L*p, rep(0,p)
      )
      res$ll
    })

    rho_ests[rep] <- rho_grid[which.max(ll_vals)]
  }

  cat(sprintf("  n_fields = %3d: mean = %.3f  sd = %.3f  bias = %.3f\n",
              n_fields, mean(rho_ests), sd(rho_ests), mean(rho_ests) - rho_true))

  results <- rbind(results, data.frame(
    n_fields = n_fields,
    n_total  = n_total,
    mean     = mean(rho_ests),
    sd       = sd(rho_ests),
    bias     = mean(rho_ests) - rho_true,
    rmse     = sqrt(mean((rho_ests - rho_true)^2))
  ))
}

# Plot: mean estimate vs n_fields
plot(results$n_fields, results$mean,
     type="b", pch=19, col="#2166ac", log="x",
     xlab="n_fields (log scale)", ylab="mean estimated rho",
     main="Mean rho estimate vs n_fields",
     ylim=c(0.4, 1.0))
abline(h=rho_true, col="red", lty=2)
segments(results$n_fields,
         results$mean - results$sd,
         results$n_fields,
         results$mean + results$sd,
         col="#2166ac", lwd=2)
legend("topright", sprintf("True rho = %.1f", rho_true),
       col="red", lty=2, bty="n")

# Plot: SD vs n_fields
plot(results$n_fields, results$sd,
     type="b", pch=19, col="#e66101", log="xy",
     xlab="n_fields (log scale)", ylab="SD (log scale)",
     main="Estimation SD vs n_fields\n(should shrink as 1/sqrt(n_fields))")
lines(results$n_fields,
      results$sd[1] * sqrt(results$n_fields[1] / results$n_fields),
      col="gray", lty=2, lwd=2)
legend("topright", c("Empirical SD", "1/sqrt(n) reference"),
       col=c("#e66101","gray"), lty=c(1,2), pch=c(19,NA), bty="n")

# Plot: bias vs n_fields
plot(results$n_fields, results$bias,
     type="b", pch=19, col="#762a83", log="x",
     xlab="n_fields (log scale)", ylab="Bias",
     main="Bias vs n_fields\n(should shrink toward 0)")
abline(h=0, col="red", lty=2)

# Plot: RMSE vs n_fields
plot(results$n_fields, results$rmse,
     type="b", pch=19, col="#4dac26", log="xy",
     xlab="n_fields (log scale)", ylab="RMSE (log scale)",
     main="RMSE vs n_fields")

par(mfrow=c(1,1))
dev.off()
cat("\nPlots saved to vignettes/multi_field_consistency.pdf\n")
