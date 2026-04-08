# Debug: why can't we identify rho?
# Key question: does AQ^{-1}(rho)A' change with rho, or is it always ~scalar*I?
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

p   <- 100
n   <- 500
set.seed(42)
A_random  <- matrix(rnorm(n*p), n, p) / sqrt(p)
A_spatial <- {
  rows <- rep(seq_len(n), each=5)
  cols <- unlist(lapply(seq_len(n), function(i) {
    start <- sample(seq_len(p-4), 1); start:(start+4)
  }))
  as.matrix(Matrix::sparseMatrix(i=rows, j=cols,
                                 x=rep(1/5,n*5), dims=c(n,p)))
}

rho_grid <- seq(0.1, 0.9, by=0.1)

pdf("vignettes/debug_identifiability.pdf", width=14, height=10)
par(mfrow=c(2,4))

# -----------------------------------------------------------------------
# 1. How much does tr(Q^{-1}(rho)) vary with rho?
# -----------------------------------------------------------------------
tr_Qinv <- sapply(rho_grid, function(rho) {
  Q <- make_SAR_Q(p, rho)
  sum(diag(solve(as.matrix(Q))))
})
plot(rho_grid, tr_Qinv, type="b", pch=19, col="#2166ac",
     xlab="rho", ylab="tr(Q^{-1})",
     main="1. tr(Q^{-1}) vs rho\n(if this varies, signal exists)")

# -----------------------------------------------------------------------
# 2. How much does A Q^{-1} A' vary with rho for random A?
# -----------------------------------------------------------------------
eigen_AQA_random <- t(sapply(rho_grid, function(rho) {
  Q    <- make_SAR_Q(p, rho)
  Qinv <- solve(as.matrix(Q))
  AQA  <- A_random %*% Qinv %*% t(A_random)
  eig  <- eigen(AQA, symmetric=TRUE, only.values=TRUE)$values
  c(min=min(eig), max=max(eig), tr=sum(eig))
}))

plot(rho_grid, eigen_AQA_random[,"tr"], type="b", pch=19, col="#e66101",
     xlab="rho", ylab="tr(AQ^{-1}A')",
     main="2. tr(AQ^{-1}A') random A\n(should vary if identifiable)")

plot(rho_grid, eigen_AQA_random[,"max"]/eigen_AQA_random[,"min"],
     type="b", pch=19, col="#e66101",
     xlab="rho", ylab="max/min eigenvalue",
     main="3. Condition of AQ^{-1}A' random A\n(1=isotropic, flat ll)")

# -----------------------------------------------------------------------
# 3. Same for spatial A
# -----------------------------------------------------------------------
eigen_AQA_spatial <- t(sapply(rho_grid, function(rho) {
  Q    <- make_SAR_Q(p, rho)
  Qinv <- solve(as.matrix(Q))
  AQA  <- A_spatial %*% Qinv %*% t(A_spatial)
  eig  <- eigen(AQA, symmetric=TRUE, only.values=TRUE)$values
  c(min=min(eig), max=max(eig), tr=sum(eig))
}))

plot(rho_grid, eigen_AQA_spatial[,"tr"], type="b", pch=19, col="#4dac26",
     xlab="rho", ylab="tr(AQ^{-1}A')",
     main="4. tr(AQ^{-1}A') spatial A\n(should vary if identifiable)")

plot(rho_grid, eigen_AQA_spatial[,"max"]/eigen_AQA_spatial[,"min"],
     type="b", pch=19, col="#4dac26",
     xlab="rho", ylab="max/min eigenvalue",
     main="5. Condition of AQ^{-1}A' spatial A")

# -----------------------------------------------------------------------
# 4. What A WOULD make rho identifiable?
# Use A = eigenvectors of Q -- perfect alignment
# -----------------------------------------------------------------------
Q_true  <- make_SAR_Q(p, 0.7)
eig_Q   <- eigen(as.matrix(Q_true), symmetric=TRUE)
A_eigen <- t(eig_Q$vectors[, 1:min(n,p)])   # n x p, rows = eigenvectors

eigen_AQA_eigen <- t(sapply(rho_grid, function(rho) {
  Q    <- make_SAR_Q(p, rho)
  Qinv <- solve(as.matrix(Q))
  AQA  <- A_eigen %*% Qinv %*% t(A_eigen)
  eig  <- eigen(AQA, symmetric=TRUE, only.values=TRUE)$values
  c(min=min(eig), max=max(eig), tr=sum(eig))
}))

plot(rho_grid, eigen_AQA_eigen[,"tr"], type="b", pch=19, col="#762a83",
     xlab="rho", ylab="tr(AQ^{-1}A')",
     main="6. tr(AQ^{-1}A') eigenvector A\n(maximum identifiability)")

# -----------------------------------------------------------------------
# 5. Compare log|Sigma_y| across designs
# Sigma_y = sigma2e*I + sigma2b * A Q^{-1} A'
# -----------------------------------------------------------------------
sigma2e <- 0.5; sigma2b <- 4
logdet_Sigmay <- function(A_mat, rho) {
  Q    <- make_SAR_Q(p, rho)
  Qinv <- solve(as.matrix(Q))
  Sy   <- sigma2e * diag(nrow(A_mat)) + sigma2b * A_mat %*% Qinv %*% t(A_mat)
  as.numeric(determinant(Sy, logarithm=TRUE)$modulus)
}

ld_random  <- sapply(rho_grid, function(r) logdet_Sigmay(A_random,  r))
ld_spatial <- sapply(rho_grid, function(r) logdet_Sigmay(A_spatial, r))

plot(rho_grid, ld_random - ld_random[1],
     type="b", pch=19, col="#e66101",
     xlab="rho", ylab="log|Sigma_y| - baseline",
     main="7. Change in log|Sigma_y| vs rho\n(flat = not identifiable)",
     ylim=range(c(ld_random-ld_random[1], ld_spatial-ld_spatial[1])))
lines(rho_grid, ld_spatial - ld_spatial[1], type="b", pch=19, col="#4dac26")
legend("topright", c("random A","spatial A"), col=c("#e66101","#4dac26"),
       pch=19, lty=1, bty="n")

par(mfrow=c(1,1))
dev.off()
cat("Plots saved to vignettes/debug_identifiability.pdf\n")
