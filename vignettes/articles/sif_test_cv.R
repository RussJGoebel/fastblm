library(downscaling)
library(Matrix)
pkgload::load_all(".")

# --- data ---------------------------------------------------------------
A <- compute_A_matrix(target_grid, soundings)
y <- soundings$SIF_757nm
W <- compute_W_matrix(target_grid, grid_shape = "other")

keep <- !is.na(y)
y    <- y[keep]
A    <- A[keep, ]
n    <- nrow(A); p <- ncol(A)
cat(sprintf("n = %d, p = %d\n", n, p))

# --- shared Q_fun -------------------------------------------------------
Q_fun <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(p) - rho * W
  Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
  ld  <- 2 * as.numeric(Matrix::determinant(S, logarithm = TRUE)$modulus)
  list(Q = Q, log_det_Q = ld, Q_matrix = Q)
}

# --- REML ---------------------------------------------------------------
cat("\n=== REML ===\n")
tuned_reml <- tune_reml(
  y             = y,
  A             = A,
  Q_fun         = Q_fun,
  theta_init    = c(rho = 0.5),
  lower         = c(rho = 0.01),
  upper         = c(rho = 0.999),
  logdet_method = "lanczos",
  verbose       = TRUE
)
print(tuned_reml)

# --- CV -----------------------------------------------------------------
cat("\n=== 10-fold CV ===\n")
tuned_cv <- tune_cv(
  y          = y,
  A          = A,
  Q_fun      = Q_fun,
  theta_init = c(rho = 0.5),
  lower      = c(rho = 0.01),
  upper      = c(rho = 0.999),
  k          = 10L,
  score      = "mse",
  verbose    = TRUE
)
print(tuned_cv)

# --- comparison ---------------------------------------------------------
cat("\n=== COMPARISON ===\n")
cat(sprintf("         %10s  %10s\n", "REML", "CV"))
cat(sprintf("rho      %10.4f  %10.4f\n", tuned_reml$theta[["rho"]], tuned_cv$theta[["rho"]]))
cat(sprintf("phi      %10.4f  %10.4f\n", tuned_reml$phi,            tuned_cv$phi))
cat(sprintf("sigma2e  %10.4f  %10.4f\n", tuned_reml$sigma2e,        tuned_cv$sigma2e))
cat(sprintf("sigma2b  %10.4f  %10.4f\n", tuned_reml$sigma2b,        tuned_cv$sigma2b))
cat(sprintf("tau      %10.4f  %10.4f\n",
            1 / (tuned_reml$phi * tuned_reml$sigma2e),
            1 / (tuned_cv$phi   * tuned_cv$sigma2e)))
