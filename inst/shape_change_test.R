# Does AQ^{-1}A' change shape with rho for different A types?
# And does the marginal ll properly exploit this?
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
rho_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# -----------------------------------------------------------------------
# Different A types
# -----------------------------------------------------------------------
set.seed(42)

# 1. Random A
A_random <- matrix(rnorm(100*p), 100, p) / sqrt(p)

# 2. Smooth local A -- each obs is weighted average of nearby basis functions
make_smooth_A <- function(n, p, bw=3/p) {
  obs  <- seq(0,1,length.out=n)
  bas  <- seq(0,1,length.out=p)
  A    <- matrix(0, n, p)
  for (i in seq_len(n)) {
    w <- exp(-0.5*((obs[i]-bas)/bw)^2)
    A[i,] <- w/sum(w)
  }
  A
}
A_smooth <- make_smooth_A(100, p)

# 3. Non-overlapping blocks -- each obs averages one contiguous block
# Like a 3x3 grid cell average -- hypothesis: NOT identifiable
A_blocks <- matrix(0, p, p)  # n=p, each obs = one non-overlapping block of 1
for (i in seq_len(p)) A_blocks[i, i] <- 1  # identity -- direct observation
# Actually try non-overlapping blocks of size 5
n_blocks <- p %/% 5
A_nonoverlap <- matrix(0, n_blocks, p)
for (i in seq_len(n_blocks)) {
  A_nonoverlap[i, ((i-1)*5+1):(i*5)] <- 1/5
}

# 4. Overlapping blocks of size 5, step 1 -- dense overlap
n_overlap <- p - 4
A_overlap <- matrix(0, n_overlap, p)
for (i in seq_len(n_overlap)) {
  A_overlap[i, i:(i+4)] <- 1/5
}

pdf("vignettes/shape_change_test.pdf", width=14, height=12)
par(mfrow=c(3,4))

A_list  <- list(random=A_random, smooth=A_smooth,
                nonoverlap=A_nonoverlap, overlap=A_overlap)
A_names <- names(A_list)

# -----------------------------------------------------------------------
# For each A: plot how AQ^{-1}A' changes with rho
# Specifically: plot the off-diagonal structure at rho_low vs rho_high
# -----------------------------------------------------------------------
for (nm in A_names) {
  A_mat <- A_list[[nm]]
  n_mat <- nrow(A_mat)

  # compute AQ^{-1}A' at two rho values
  Q_low  <- make_SAR_Q(p, 0.1)
  Q_high <- make_SAR_Q(p, 0.9)
  AQA_low  <- A_mat %*% solve(as.matrix(Q_low))  %*% t(A_mat)
  AQA_high <- A_mat %*% solve(as.matrix(Q_high)) %*% t(A_mat)

  # how much does it change?
  change <- AQA_high - AQA_low
  rel_change <- norm(change, "F") / norm(AQA_low, "F")

  # plot first row of AQA to see covariance structure
  plot(seq_len(n_mat), AQA_low[1,],
       type="l", col="#2166ac", lwd=2,
       xlab="observation j", ylab="Cov(y_1, y_j)",
       main=sprintf("%s\nCov structure (row 1 of AQ^{-1}A')", nm),
       ylim=range(c(AQA_low[1,], AQA_high[1,])))
  lines(seq_len(n_mat), AQA_high[1,], col="#e66101", lwd=2)
  legend("topright", c("rho=0.1","rho=0.9"),
         col=c("#2166ac","#e66101"), lty=1, lwd=2, bty="n", cex=0.8)
  title(sub=sprintf("rel change in AQA': %.3f", rel_change), cex.sub=0.8)
}

# -----------------------------------------------------------------------
# For each A: plot marginal ll surface over rho
# -----------------------------------------------------------------------
phi_true    <- 5
sigma2e_true <- 0.5
Q_true <- make_SAR_Q(p, rho_true)
CQ     <- Matrix::Cholesky(Q_true)
probes <- matrix(sample(c(-1L,1L), p*50, replace=TRUE), p, 50)

for (nm in A_names) {
  A_mat <- A_list[[nm]]
  n_mat <- nrow(A_mat)
  set.seed(42)
  x_true <- as.numeric(Matrix::solve(CQ, rnorm(p))) * sqrt(phi_true*sigma2e_true)
  y      <- as.numeric(A_mat %*% x_true + rnorm(n_mat, sd=sqrt(sigma2e_true)))

  # exact marginal ll
  rho_fine <- seq(0.05, 0.90, by=0.05)
  ll_vals  <- sapply(rho_fine, function(rho) {
    Q    <- make_SAR_Q(p, rho)
    Qinv <- solve(as.matrix(Q))
    Sy   <- sigma2e_true * (diag(n_mat) + phi_true * A_mat %*% Qinv %*% t(A_mat))
    ld   <- as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
    qf   <- as.numeric(t(y) %*% solve(Sy, y))
    -n_mat/2*log(2*pi) - 1/2*ld - 1/2*qf
  })

  # also profile phi
  ll_profiled <- sapply(rho_fine, function(rho) {
    Q    <- make_SAR_Q(p, rho)
    Qinv <- solve(as.matrix(Q))
    Sy_shape <- diag(n_mat) + phi_true * A_mat %*% Qinv %*% t(A_mat)

    # try a range of phi
    best <- -Inf
    for (phi in seq(0.5, 20, by=0.5)) {
      Sy <- (diag(n_mat) + phi * A_mat %*% Qinv %*% t(A_mat))
      s2 <- as.numeric(t(y) %*% solve(Sy, y)) / n_mat
      if (s2 <= 0) next
      ld <- as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
      ll <- -n_mat/2*log(s2) - 1/2*ld - n_mat/2
      if (ll > best) best <- ll
    }
    best
  })

  peak_fixed  <- rho_fine[which.max(ll_vals)]
  peak_profiled <- rho_fine[which.max(ll_profiled)]

  plot(rho_fine, ll_profiled - max(ll_profiled),
       type="b", pch=19, col="#4dac26",
       xlab="rho", ylab="ll (normalized)",
       main=sprintf("%s\nMarginal ll (phi profiled)", nm))
  lines(rho_fine, ll_vals - max(ll_vals),
        type="b", pch=19, col="#2166ac", cex=0.5)
  abline(v=rho_true,      col="red",     lty=2)
  abline(v=peak_profiled, col="#4dac26", lty=2)
  legend("bottomleft",
         c(sprintf("profiled phi (peak=%.2f)", peak_profiled),
           sprintf("fixed true phi (peak=%.2f)", peak_fixed),
           sprintf("true rho=%.1f", rho_true)),
         col=c("#4dac26","#2166ac","red"),
         lty=c(1,1,2), pch=c(19,19,NA), bty="n", cex=0.7)
}

# -----------------------------------------------------------------------
# Summary: how much does AQA' change shape vs scale?
# Decompose change into trace (scale) vs off-diagonal (shape)
# -----------------------------------------------------------------------
cat("\n=== Shape vs scale change in AQ^{-1}A' ===\n")
cat(sprintf("%-12s | %-10s %-10s %-10s\n", "A type", "tr change", "offdiag", "rel_shape"))
cat(strrep("-", 50), "\n")

for (nm in A_names) {
  A_mat   <- A_list[[nm]]
  Q_low   <- make_SAR_Q(p, 0.1)
  Q_high  <- make_SAR_Q(p, 0.9)
  AQA_low  <- A_mat %*% solve(as.matrix(Q_low))  %*% t(A_mat)
  AQA_high <- A_mat %*% solve(as.matrix(Q_high)) %*% t(A_mat)

  change      <- AQA_high - AQA_low
  tr_change   <- sum(diag(change)) / sum(diag(AQA_low))
  # off diagonal change relative to diagonal change
  diag_change <- sum(diag(change)^2)
  total_change <- sum(change^2)
  shape_change <- sqrt((total_change - diag_change) / total_change)

  cat(sprintf("%-12s | %-10.3f %-10.3f %-10.3f\n",
              nm, tr_change, sqrt(sum(change[row(change)!=col(change)]^2)),
              shape_change))
}

par(mfrow=c(1,1))
dev.off()
cat("Plots saved to vignettes/shape_change_test.pdf\n")
