# =============================================================================
# sif_penalty_test.R
#
# Fits the SIF downscaling model with and without a soft RSR penalty,
# using both REML and CV. Compares rho, phi/tau, and lambda estimates.
#
# The soft penalty adds lambda * t(C) %*% C to the SAR precision, where
# C = t(X_obs) %*% A is the RSR constraint matrix. As lambda -> Inf this
# recovers the hard RSR constraint.
#
# Usage:
#   library(downscaling); pkgload::load_all(".")
#   source("utils.R"); source("solvers.R"); source("fit.R")
#   source("reml.R"); source("cv.R")
#   source("sif_penalty_test.R")
# =============================================================================

library(downscaling)
library(Matrix)
pkgload::load_all(".")

# --- data -------------------------------------------------------------------
A <- compute_A_matrix(target_grid, soundings)
y <- soundings$SIF_757nm
W <- compute_W_matrix(target_grid, grid_shape = "other")

keep <- !is.na(y)
y    <- y[keep]
A    <- A[keep, ]
n    <- nrow(A)
p    <- ncol(A)
cat(sprintf("n = %d, p = %d\n", n, p))

# --- covariates and constraint ----------------------------------------------
# Subset covariates by keep before cbind to ensure conformable dimensions
X_obs <- cbind(1, soundings$water_fraction)[keep, , drop = FALSE]
X_obs <- as.matrix(X_obs)
q     <- ncol(X_obs)
C_mat <- as.matrix(t(X_obs) %*% A)    # q x p
CtC   <- Matrix::crossprod(C_mat)      # p x p, precompute once

cat(sprintf("Covariates: %d x %d  (intercept + water fraction)\n", nrow(X_obs), q))
cat(sprintf("Constraint C: %d x %d\n\n", nrow(C_mat), ncol(C_mat)))

# --- Q_fun factories --------------------------------------------------------

# Base: SAR only, theta = c(rho)
make_Q_fun_base <- function() {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(p) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    ld  <- 2 * as.numeric(Matrix::determinant(S, logarithm = TRUE)$modulus)
    list(Q = Q, log_det_Q = ld, Q_matrix = Q)
  }
}

# Penalised: SAR + lambda * CtC, theta = c(rho, log_lambda)
make_Q_fun_penalised <- function() {
  function(theta) {
    rho <- theta[["rho"]]
    lam <- exp(theta[["log_lambda"]])
    S   <- Matrix::Diagonal(p) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S) + lam * CtC)
    ld  <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
    list(Q = Q, log_det_Q = ld, Q_matrix = Q)
  }
}

Q_fun_base <- make_Q_fun_base()
Q_fun_pen  <- make_Q_fun_penalised()

# --- REML: base -------------------------------------------------------------
cat("=== REML: rho only ===\n")
tuned_reml_base <- tune_reml(
  y             = y,
  A             = A,
  Q_fun         = Q_fun_base,
  theta_init    = c(rho = 0.5),
  lower         = c(rho = 0.01),
  upper         = c(rho = 0.999),
  logdet_method = "lanczos",
  verbose       = TRUE
)
print(tuned_reml_base)

# --- REML: rho + lambda -----------------------------------------------------
cat("\n=== REML: rho + lambda ===\n")
tuned_reml_pen <- tune_reml(
  y             = y,
  A             = A,
  Q_fun         = Q_fun_pen,
  theta_init    = c(rho = 0.5, log_lambda = log(100)),
  lower         = c(rho = 0.01, log_lambda = log(1)),
  upper         = c(rho = 0.999, log_lambda = log(1e8)),
  logdet_method = "lanczos",
  verbose       = TRUE
)
print(tuned_reml_pen)

# --- CV: base ---------------------------------------------------------------
cat("\n=== CV: rho only ===\n")
tuned_cv_base <- tune_cv(
  y             = y,
  A             = A,
  Q_fun         = Q_fun_base,
  theta_init    = c(rho = 0.5),
  lower         = c(rho = 0.01),
  upper         = c(rho = 0.999),
  k             = 5L,
  score         = "mse",
  log_phi_upper = log(10),
  seed          = 1L,
  verbose       = TRUE
)
print(tuned_cv_base)

# --- CV: rho + lambda -------------------------------------------------------
cat("\n=== CV: rho + lambda ===\n")
tuned_cv_pen <- tune_cv(
  y             = y,
  A             = A,
  Q_fun         = Q_fun_pen,
  theta_init    = c(rho = 0.5, log_lambda = log(100)),
  lower         = c(rho = 0.01, log_lambda = log(1)),
  upper         = c(rho = 0.999, log_lambda = log(1e8)),
  k             = 5L,
  score         = "mse",
  log_phi_upper = log(10),
  seed          = 1L,
  verbose       = TRUE
)
print(tuned_cv_pen)

# --- comparison -------------------------------------------------------------
cat("\n=====================================================\n")
cat(" COMPARISON\n")
cat("=====================================================\n")

fmt <- function(label, reml_base, reml_pen, cv_base, cv_pen) {
  cat(sprintf("%-14s  %10.4f  %10.4f  %10.4f  %10.4f\n",
              label, reml_base, reml_pen, cv_base, cv_pen))
}

cat(sprintf("%-14s  %10s  %10s  %10s  %10s\n",
            "", "REML_base", "REML_pen", "CV_base", "CV_pen"))
cat(sprintf("%s\n", strrep("-", 58)))

fmt("rho",
    tuned_reml_base$theta[["rho"]],
    tuned_reml_pen$theta[["rho"]],
    tuned_cv_base$theta[["rho"]],
    tuned_cv_pen$theta[["rho"]])

fmt("lambda",
    NA,
    exp(tuned_reml_pen$theta[["log_lambda"]]),
    NA,
    exp(tuned_cv_pen$theta[["log_lambda"]]))

fmt("phi",
    tuned_reml_base$phi,
    tuned_reml_pen$phi,
    tuned_cv_base$phi,
    tuned_cv_pen$phi)

fmt("tau (= 1/phi)",
    1 / tuned_reml_base$phi,
    1 / tuned_reml_pen$phi,
    1 / tuned_cv_base$phi,
    1 / tuned_cv_pen$phi)

fmt("sigma2e",
    tuned_reml_base$sigma2e,
    tuned_reml_pen$sigma2e,
    tuned_cv_base$sigma2e,
    tuned_cv_pen$sigma2e)

cat("\nNote: lambda -> upper bound (1e8) means data prefer hard constraint.\n")
cat("      lambda in middle means soft constraint preferred.\n")
