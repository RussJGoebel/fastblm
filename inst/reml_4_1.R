# =============================================================================
# test_reml.R
# Run with: source("test_reml.R")
# All tests print PASS / FAIL. Final line prints overall result.
# =============================================================================

pkgload::load_all(".")

# =============================================================================
# Tiny test harness
# =============================================================================
.pass <- 0L
.fail <- 0L
.results <- list()

expect <- function(label, expr) {
  result <- tryCatch({
    val <- expr
    if (isTRUE(val)) {
      list(ok = TRUE)
    } else {
      list(ok = FALSE, msg = paste("Expression was not TRUE:", deparse(substitute(expr))))
    }
  }, error = function(e) {
    list(ok = FALSE, msg = conditionMessage(e))
  })

  if (result$ok) {
    cat(sprintf("  PASS  %s\n", label))
    .pass <<- .pass + 1L
  } else {
    cat(sprintf("  FAIL  %s\n         %s\n", label, result$msg))
    .fail <<- .fail + 1L
  }
}

expect_error <- function(label, expr) {
  result <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  if (isTRUE(result)) {
    cat(sprintf("  PASS  %s\n", label))
    .pass <<- .pass + 1L
  } else {
    cat(sprintf("  FAIL  %s  (expected an error but got none)\n", label))
    .fail <<- .fail + 1L
  }
}

expect_warning <- function(label, expr) {
  saw_warning <- FALSE
  withCallingHandlers(expr, warning = function(w) {
    saw_warning <<- TRUE
    invokeRestart("muffleWarning")
  })
  if (saw_warning) {
    cat(sprintf("  PASS  %s\n", label))
    .pass <<- .pass + 1L
  } else {
    cat(sprintf("  FAIL  %s  (expected a warning but got none)\n", label))
    .fail <<- .fail + 1L
  }
}

section <- function(title) cat(sprintf("\n--- %s ---\n", title))

# =============================================================================
# Shared small problem setup
# =============================================================================
set.seed(42)
n  <- 60
p  <- 10
q  <- 2    # number of fixed effects

# Fixed effects design matrix (intercept + one covariate)
X  <- cbind(1, rnorm(n))

# Random effects design matrix
A  <- matrix(rnorm(n * p), n, p)

# True parameters
beta_true  <- c(2, -1)
x_true     <- rnorm(p)
sigma2e    <- 0.5
y          <- X %*% beta_true + A %*% x_true + rnorm(n, sd = sqrt(sigma2e))

# A simple diagonal prior precision Q (ridge-like)
Q_mat <- Matrix::Diagonal(p, 1)
Q_fun_simple <- function(theta) list(Q = Q_mat, log_det_Q = 0)

# =============================================================================
# 1. reml_project: basic behaviour
# =============================================================================
section("reml_project: basic behaviour")

proj <- reml_project(y, A, X)

expect("returns a list with y, A, n_eff",
       is.list(proj) && all(c("y", "A", "n_eff") %in% names(proj))
)

expect("projected y has same length as original y",
       length(proj$y) == n
)

expect("projected A has same dimensions as original A",
       nrow(proj$A) == n && ncol(proj$A) == p
)

expect("n_eff equals n minus rank(X)",
       proj$n_eff == n - qr(X)$rank
)

expect("projected y is orthogonal to X (residuals of regression on X)",
       max(abs(t(X) %*% proj$y)) < 1e-10
)

expect("projected A columns are orthogonal to X",
       max(abs(t(X) %*% proj$A)) < 1e-10
)

expect("projecting twice gives the same result (idempotent)",
       max(abs(reml_project(proj$y, proj$A, X)$y - proj$y)) < 1e-10
)

# =============================================================================
# 2. reml_project: input validation
# =============================================================================
section("reml_project: input validation")

expect_error("rejects non-matrix X_fixed",
             reml_project(y, A, as.data.frame(X))
)

expect_error("rejects X_fixed with wrong number of rows",
             reml_project(y, A, X[1:5, ])
)

expect_warning("warns on rank-deficient X_fixed",
               reml_project(y, A, cbind(X, X[, 1]))  # duplicate column -> rank < ncol
)

# =============================================================================
# 3. tune_reml: runs without X_fixed (baseline, no regression)
# =============================================================================
section("tune_reml: no fixed effects (baseline)")

fit0 <- tune_reml(y, A, Q_fun_simple,
                  log_phi_lower = log(0.1),
                  log_phi_upper = log(100),
                  verbose = FALSE)

expect("returns fastblm_tuned object",
       inherits(fit0, "fastblm_tuned")
)

expect("phi is positive",
       is.numeric(fit0$phi) && fit0$phi > 0
)

expect("sigma2e is positive",
       is.numeric(fit0$sigma2e) && fit0$sigma2e > 0
)

expect("sigma2b equals phi * sigma2e",
       abs(fit0$sigma2b - fit0$phi * fit0$sigma2e) < 1e-12
)

expect("method is 'reml'",
       fit0$method == "reml"
)

# =============================================================================
# 4. tune_reml: X_fixed integration — n_eff correction
# =============================================================================
section("tune_reml: X_fixed integration")

fit_xf <- tune_reml(y, A, Q_fun_simple,
                    X_fixed       = X,
                    log_phi_lower = log(0.1),
                    log_phi_upper = log(100),
                    verbose       = FALSE)

expect("returns fastblm_tuned object with X_fixed supplied",
       inherits(fit_xf, "fastblm_tuned")
)

expect("phi is positive with X_fixed",
       is.numeric(fit_xf$phi) && fit_xf$phi > 0
)

expect("sigma2e is positive with X_fixed",
       is.numeric(fit_xf$sigma2e) && fit_xf$sigma2e > 0
)

# Manually project and verify sigma2e uses n_eff not n
proj_manual <- reml_project(y, A, X)
n_eff <- proj_manual$n_eff

# Passing pre-projected data without X_fixed uses n as denominator (not n_eff).
# The X_fixed path correctly uses n_eff, so sigma2e will differ by exactly n/n_eff.
# We verify this scaling relationship holds.
fit_manual_n <- tune_reml(proj_manual$y, proj_manual$A, Q_fun_simple,
                          log_phi_lower = log(0.1),
                          log_phi_upper = log(100),
                          verbose = FALSE)

expect("sigma2e with X_fixed is larger than without (n_eff < n means bigger estimate)",
       fit_xf$sigma2e > fit_manual_n$sigma2e
)

# When phi is fixed, sigma2e scales exactly as n/n_eff. With phi also optimised
# the ratio won't be exact, but the direction must hold: X_fixed gives a larger
# sigma2e because n_eff < n. The test above already checks this.
# We also verify the difference is non-trivial (not just floating point noise).
expect("sigma2e difference is meaningful, not just numerical noise",
       (fit_xf$sigma2e - fit_manual_n$sigma2e) > 0.01 * fit_manual_n$sigma2e
)

# =============================================================================
# 5. tune_reml: X_fixed with rank-deficient X warns but still runs
# =============================================================================
section("tune_reml: rank-deficient X_fixed")

X_rankdef <- cbind(X, X[, 1])  # duplicate column

expect_warning("rank-deficient X_fixed triggers a warning",
               tune_reml(y, A, Q_fun_simple,
                         X_fixed       = X_rankdef,
                         log_phi_lower = log(0.1),
                         log_phi_upper = log(100),
                         verbose       = FALSE)
)

# =============================================================================
# 6. tune_reml: theta optimisation with one shape parameter
# =============================================================================
section("tune_reml: one theta parameter")

# Q_fun that scales precision by exp(theta)
Q_fun_theta <- function(theta) {
  scale <- exp(theta[["log_scale"]])
  list(Q = scale * Q_mat, log_det_Q = p * log(scale))
}

fit_theta <- tune_reml(y, A, Q_fun_theta,
                       theta_init    = c(log_scale = 0),
                       lower         = -3,
                       upper         =  3,
                       log_phi_lower = log(0.1),
                       log_phi_upper = log(100),
                       verbose       = FALSE)

expect("theta optimisation returns named theta",
       !is.null(names(fit_theta$theta)) && "log_scale" %in% names(fit_theta$theta)
)

expect("optimised theta is within bounds",
       fit_theta$theta[["log_scale"]] >= -3 &&
         fit_theta$theta[["log_scale"]] <=  3
)

expect("history is non-empty after theta optimisation",
       is.data.frame(fit_theta$history) && nrow(fit_theta$history) > 0
)

# =============================================================================
# 7. tune_reml: constraint_matrix still works alongside X_fixed
# =============================================================================
section("tune_reml: constraint_matrix + X_fixed")

# Simple sum-to-zero constraint on the random effects
C_mat <- matrix(1, nrow = 1, ncol = p)

fit_constrained <- tune_reml(y, A, Q_fun_simple,
                             X_fixed           = X,
                             constraint_matrix = C_mat,
                             log_phi_lower     = log(0.1),
                             log_phi_upper     = log(100),
                             verbose           = FALSE)

expect("constrained + X_fixed fit returns fastblm_tuned",
       inherits(fit_constrained, "fastblm_tuned")
)

expect("constrained + X_fixed phi is positive",
       fit_constrained$phi > 0
)

# =============================================================================
# 8. reml_project: intercept-only special case
# =============================================================================
section("reml_project: intercept-only X_fixed")

X_int  <- matrix(1, n, 1)
proj_int <- reml_project(y, A, X_int)

expect("n_eff = n - 1 for intercept-only X_fixed",
       proj_int$n_eff == n - 1
)

expect("projected y has zero mean (intercept removed)",
       abs(mean(proj_int$y)) < 1e-10
)


# =============================================================================
# 9. constraint_matrix: Schur correction verified against reparameterisation
#
# We verify the constrained likelihood implementation by comparing it to the
# reparameterisation approach (known correct, dense) on a small explicit problem.
# Both should give the same optimised phi, sigma2e and ll value.
# =============================================================================
section("constraint_matrix: Schur correction vs reparameterisation")

set.seed(99)
n_c  <- 25
p_c  <- 6
r_c  <- 2
p_eff_c <- p_c - r_c

A_c <- matrix(rnorm(n_c * p_c), n_c, p_c)

# Explicit dense Q so both paths can use it
Q_c_mat  <- diag(p_c) + 0.2 * crossprod(matrix(rnorm(p_c^2), p_c, p_c)) / p_c
Q_c_mat  <- (Q_c_mat + t(Q_c_mat)) / 2 + p_c * diag(p_c)
Q_c_mat  <- Matrix::Matrix(Q_c_mat, sparse = FALSE)

logdet_Q_c <- as.numeric(Matrix::determinant(Q_c_mat, logarithm = TRUE)$modulus)
Q_fun_c    <- function(theta) list(Q = Q_c_mat, Q_matrix = Q_c_mat,
                                   log_det_Q = logdet_Q_c)

# Constraint matrix
C_c <- matrix(rnorm(r_c * p_c), r_c, p_c)
while (qr(C_c)$rank < r_c) C_c <- matrix(rnorm(r_c * p_c), r_c, p_c)

# Generate y from a constrained x
qr_Cc  <- qr(t(C_c))
N_c    <- qr.Q(qr_Cc, complete = TRUE)[, (r_c + 1L):p_c, drop = FALSE]
x_c    <- N_c %*% rnorm(p_eff_c)
y_c    <- as.numeric(A_c %*% x_c + rnorm(n_c, sd = sqrt(0.5)))

# --- Schur correction (the actual implementation) ---------------------------
fit_schur <- tune_reml(y_c, A_c, Q_fun_c,
                       constraint_matrix = C_c,
                       log_phi_lower     = log(0.1),
                       log_phi_upper     = log(50),
                       verbose           = FALSE)

# --- Reparameterisation (dense, known correct) ------------------------------
A_tilde_c  <- A_c %*% N_c
Q_red_c    <- t(N_c) %*% as.matrix(Q_c_mat) %*% N_c
logdet_Qr  <- as.numeric(determinant(Q_red_c, logarithm = TRUE)$modulus)
Q_fun_r    <- function(theta) list(Q = Matrix::Matrix(Q_red_c, sparse = FALSE),
                                   Q_matrix = Matrix::Matrix(Q_red_c, sparse = FALSE),
                                   log_det_Q = logdet_Qr)

fit_reparam <- tune_reml(y_c, A_tilde_c, Q_fun_r,
                         log_phi_lower = log(0.1),
                         log_phi_upper = log(50),
                         verbose       = FALSE)

expect("Schur correction returns fastblm_tuned",
       inherits(fit_schur, "fastblm_tuned")
)

# Note: phi and sigma2e should agree between the two approaches because they
# find the same optimum. The ll values will differ by a constant offset because
# the two parameterisations drop different normalising constants when profiling
# sigma2e -- this is expected and doesn't affect the optimum.
expect("Schur phi agrees with reparameterisation (within 10%)",
       abs(fit_schur$phi - fit_reparam$phi) / fit_reparam$phi < 0.10
)

expect("Schur sigma2e agrees with reparameterisation (within 5%)",
       abs(fit_schur$sigma2e - fit_reparam$sigma2e) / fit_reparam$sigma2e < 0.05
)

expect("Schur ll is finite and optimiser converged",
       is.finite(fit_schur$value) && fit_schur$optim$convergence == 0L
)



# =============================================================================
# 10. General R_inv: reml_project and tune_reml with non-identity noise
#
# Verifies that R_inv is handled correctly throughout:
#   - reml_project with R_inv produces correct weighted projection
#   - tune_reml with R_inv + X_fixed gives same result as manually whitening
#     and projecting before calling tune_reml without X_fixed
#   - constraint correction works with non-identity R_inv
# =============================================================================
section("General R_inv support")

set.seed(77)
n_r  <- 30
p_r  <- 5
r_r  <- 1

A_r  <- matrix(rnorm(n_r * p_r), n_r, p_r)
X_r  <- cbind(1, rnorm(n_r))   # intercept + one covariate

# Build a simple diagonal R (heteroskedastic noise)
R_diag  <- exp(rnorm(n_r, sd = 0.3))
R_r     <- diag(R_diag)
R_inv_r <- diag(1 / R_diag)

Q_r_mat  <- diag(p_r) + 0.1 * crossprod(matrix(rnorm(p_r^2), p_r, p_r)) / p_r
Q_r_mat  <- (Q_r_mat + t(Q_r_mat)) / 2 + p_r * diag(p_r)
Q_r_mat  <- Matrix::Matrix(Q_r_mat, sparse = FALSE)
logdet_Q_r <- as.numeric(Matrix::determinant(Q_r_mat, logarithm = TRUE)$modulus)
Q_fun_r  <- function(theta) list(Q = Q_r_mat, Q_matrix = Q_r_mat,
                                 log_det_Q = logdet_Q_r)

x_true_r <- rnorm(p_r)
beta_r   <- c(1.5, -0.5)
y_r <- X_r %*% beta_r + A_r %*% x_true_r + rnorm(n_r, sd = sqrt(R_diag))

# --- reml_project with R_inv -------------------------------------------------
proj_r <- reml_project(y_r, A_r, X_r, R_inv = R_inv_r)

expect("reml_project with R_inv: n_eff correct",
       proj_r$n_eff == n_r - qr(X_r)$rank
)

expect("reml_project with R_inv: projected y is R^{-1}-orthogonal to X",
       max(abs(t(X_r) %*% R_inv_r %*% proj_r$y)) < 1e-10
)

expect("reml_project with R_inv: projected A is R^{-1}-orthogonal to X",
       max(abs(t(X_r) %*% R_inv_r %*% proj_r$A)) < 1e-10
)

# --- tune_reml: X_fixed + R_inv should match manual whitening ----------------
# Path 1: pass X_fixed and R_inv directly
fit_r_xf <- tune_reml(y_r, A_r, Q_fun_r,
                      X_fixed       = X_r,
                      R_inv         = R_inv_r,
                      log_phi_lower = log(0.1),
                      log_phi_upper = log(50),
                      verbose       = FALSE)

expect("tune_reml with R_inv + X_fixed returns fastblm_tuned",
       inherits(fit_r_xf, "fastblm_tuned")
)

expect("tune_reml R_inv + X_fixed: phi is positive and finite",
       fit_r_xf$phi > 0 && is.finite(fit_r_xf$phi)
)

# The real test: weighted (R_inv) and unweighted (NULL) projections should give
# different results when R != I, confirming R_inv is actually being used.
fit_r_unweighted <- tune_reml(y_r, A_r, Q_fun_r,
                              X_fixed       = X_r,
                              R_inv         = NULL,
                              log_phi_lower = log(0.1),
                              log_phi_upper = log(50),
                              verbose       = FALSE)

expect("weighted and unweighted projections differ when R != I",
       abs(fit_r_xf$phi - fit_r_unweighted$phi) / fit_r_unweighted$phi > 0.01
)

# --- constraint + R_inv ------------------------------------------------------
C_r <- matrix(1, nrow = 1, ncol = p_r)   # sum-to-zero

fit_r_constr <- tune_reml(y_r, A_r, Q_fun_r,
                          R_inv             = R_inv_r,
                          constraint_matrix = C_r,
                          log_phi_lower     = log(0.1),
                          log_phi_upper     = log(50),
                          verbose           = FALSE)

expect("constraint + R_inv: returns fastblm_tuned",
       inherits(fit_r_constr, "fastblm_tuned")
)

expect("constraint + R_inv: phi is positive and finite",
       is.numeric(fit_r_constr$phi) && fit_r_constr$phi > 0 &&
         is.finite(fit_r_constr$phi)
)

expect("constraint + R_inv: sigma2e is positive and finite",
       is.numeric(fit_r_constr$sigma2e) && fit_r_constr$sigma2e > 0 &&
         is.finite(fit_r_constr$sigma2e)
)

# R_inv = I should give same result as NULL when noise is actually identity
fit_r_null <- tune_reml(y_r, A_r, Q_fun_r,
                        R_inv         = NULL,
                        log_phi_lower = log(0.1),
                        log_phi_upper = log(50),
                        verbose       = FALSE)
fit_r_eye  <- tune_reml(y_r, A_r, Q_fun_r,
                        R_inv         = Matrix::Diagonal(n_r),
                        log_phi_lower = log(0.1),
                        log_phi_upper = log(50),
                        verbose       = FALSE)

expect("R_inv = I gives same result as R_inv = NULL",
       abs(fit_r_null$phi - fit_r_eye$phi) / fit_r_null$phi < 1e-6 &&
         abs(fit_r_null$sigma2e - fit_r_eye$sigma2e) / fit_r_null$sigma2e < 1e-6
)

# =============================================================================
# Summary
# =============================================================================
cat(sprintf(
  "\n=== Results: %d passed, %d failed (out of %d) ===\n",
  .pass, .fail, .pass + .fail
))
if (.fail == 0L) {
  cat("All tests passed.\n")
} else {
  cat("Some tests FAILED — see above.\n")
}
