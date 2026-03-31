library(downscaling)
library(Matrix)
pkgload::load_all(".")

# --- data ---------------------------------------------------------------
A    <- compute_A_matrix(target_grid, soundings)
y    <- soundings$SIF_757nm
W    <- compute_W_matrix(target_grid, grid_shape = "other")

keep <- !is.na(y)
y    <- y[keep]
A    <- A[keep, ]
n    <- nrow(A); p <- ncol(A)
cat(sprintf("n = %d, p = %d\n", n, p))

# --- Q_fun with exact log|Q| -------------------------------------------
Q_fun <- function(theta) {
  rho <- theta[["rho"]]
  S   <- Matrix::Diagonal(p) - rho * W
  Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
  ld  <- 2 * as.numeric(Matrix::determinant(S, logarithm = TRUE)$modulus)
  list(Q = Q, log_det_Q = ld, Q_matrix = Q)
}

# --- tune ---------------------------------------------------------------
tuned <- tune_reml(
  y             = y,
  A             = A,
  Q_fun         = Q_fun,
  theta_init    = c(rho = 0.5),
  lower         = c(rho = 0.01),
  upper         = c(rho = 0.9999),
  logdet_method = "lanczos",   # for logdet_K only; logdet_Q is exact above
  verbose       = TRUE
)
print(tuned)

# --- fit ----------------------------------------------------------------
fit <- fit_fastblm(y, A, tuned$Q, phi = tuned$phi, solver = "cholesky")
print(fit)
