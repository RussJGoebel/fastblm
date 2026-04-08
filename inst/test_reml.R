# =============================================================================
# test_reml.R
# Run with: source("test_reml.R")
#
# This test file targets the exact-REML implementation in reml.R.
# It does NOT assume that REML is implemented by Euclidean projection of y and A.
# Instead, it tests the exact marginal objective directly.
# =============================================================================

pkgload::load_all(".")

# =============================================================================
# Tiny test harness
# =============================================================================
.pass <- 0L
.fail <- 0L

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

section <- function(title) cat(sprintf("\n--- %s ---\n", title))

# =============================================================================
# Dense reference REML helper
# =============================================================================
#
# This helper computes the dense restricted log-likelihood for the model
#
#   y = A x + X beta + e
#   e ~ N(0, sigma2e R)
#   x ~ N(0, sigma2e phi Q^{-1})
#
# optionally with the exact linear constraint C x = 0.
#
# The value returned as ll_drop_R_const matches the objective used by reml.R,
# which drops the additive constant -(1/2) log|R| when R is fixed.
# =============================================================================
dense_reml_ll <- function(y, A, Q, phi, X_fixed = NULL, R = NULL, C = NULL) {
  y <- as.numeric(y)
  A <- as.matrix(A)
  Q <- as.matrix(Q)

  n <- length(y)

  if (is.null(R)) {
    R <- diag(n)
  } else {
    R <- as.matrix(R)
  }

  Q_inv <- solve(Q)

  if (is.null(C)) {
    S <- Q_inv
  } else {
    C <- as.matrix(C)
    CQinvCt <- C %*% Q_inv %*% t(C)
    S <- Q_inv - Q_inv %*% t(C) %*% solve(CQinvCt) %*% C %*% Q_inv
  }

  H <- R + phi * A %*% S %*% t(A)
  H <- 0.5 * (H + t(H))
  H_inv <- solve(H)

  if (is.null(X_fixed)) {
    n_eff <- n
    logdet_XHX <- 0
    Hinv_y <- H_inv %*% y
    quad <- drop(crossprod(y, Hinv_y))
  } else {
    X_fixed <- as.matrix(X_fixed)

    Hinv_y <- H_inv %*% y
    Hinv_X <- H_inv %*% X_fixed

    XtHinvX <- crossprod(X_fixed, Hinv_X)
    XtHinvX <- 0.5 * (XtHinvX + t(XtHinvX))

    b <- as.numeric(crossprod(X_fixed, Hinv_y))
    alpha <- solve(XtHinvX, b)

    quad <- drop(crossprod(y, Hinv_y)) - drop(crossprod(b, alpha))
    n_eff <- n - qr(X_fixed)$rank
    logdet_XHX <- as.numeric(determinant(XtHinvX, logarithm = TRUE)$modulus)
  }

  logdet_H <- as.numeric(determinant(H, logarithm = TRUE)$modulus)
  logdet_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)

  ll_full <- -0.5 * logdet_H -
    0.5 * logdet_XHX -
    (n_eff / 2) * log(quad / n_eff)

  ll_drop_R_const <- ll_full + 0.5 * logdet_R

  list(
    ll_full = ll_full,
    ll_drop_R_const = ll_drop_R_const,
    quad = quad,
    n_eff = n_eff,
    logdet_R = logdet_R
  )
}

# =============================================================================
# Internal helper wrapper for fixed-phi evaluation
# =============================================================================
eval_new_ll <- function(y, A, Q_fun, phi, X_fixed = NULL, R_inv = NULL, C = NULL) {
  y <- as.numeric(y)
  A <- as.matrix(A)

  n <- length(y)
  p <- ncol(A)

  Rinv_obj <- .make_Rinv(R_inv, n)
  apply_Rinv <- Rinv_obj$apply

  if (!is.null(Rinv_obj$matrix)) {
    AtRinvA_matrix <- crossprod(A, Rinv_obj$matrix %*% A)
  } else {
    RinvA <- apply(A, 2, apply_Rinv)
    RinvA <- as.matrix(RinvA)
    AtRinvA_matrix <- crossprod(A, RinvA)
  }

  AtRinvA_apply <- function(v) as.numeric(AtRinvA_matrix %*% v)

  prior <- Q_fun(numeric(0))
  apply_Q <- .as_apply(prior$Q)

  set.seed(1)
  probes <- matrix(sample(c(-1L, 1L), p * 10, replace = TRUE), nrow = p)

  .eval_reml_ll(
    y = y,
    A = A,
    X_fixed = X_fixed,
    phi = phi,
    prior = prior,
    apply_Q = apply_Q,
    logdet_Q = prior$log_det_Q,
    apply_Rinv = apply_Rinv,
    AtRinvA_apply = AtRinvA_apply,
    AtRinvA_matrix = AtRinvA_matrix,
    probes = probes,
    n_lanczos_steps = 10L,
    logdet_method = "cholesky",
    pcg_tol = 1e-10,
    pcg_maxit = 1000L,
    constraint_matrix = C
  )
}

# =============================================================================
# Shared small problem setup
# =============================================================================
set.seed(42)
n  <- 60
p  <- 10
q  <- 2

X  <- cbind(1, rnorm(n))
A  <- matrix(rnorm(n * p), n, p)

beta_true <- c(2, -1)
x_true    <- rnorm(p)
sigma2e_true <- 0.5

y <- as.numeric(X %*% beta_true + A %*% x_true + rnorm(n, sd = sqrt(sigma2e_true)))

Q_mat <- Matrix::Diagonal(p, 1)
Q_fun_simple <- function(theta) {
  list(
    Q = Q_mat,
    Q_matrix = Q_mat,
    log_det_Q = 0
  )
}

# =============================================================================
# 1. Baseline fit without fixed effects
# =============================================================================
section("tune_reml: no fixed effects")

fit0 <- tune_reml(
  y, A, Q_fun_simple,
  log_phi_lower = log(0.1),
  log_phi_upper = log(100),
  verbose = FALSE
)

expect("returns fastblm_tuned object",
       inherits(fit0, "fastblm_tuned")
)

expect("phi is positive",
       is.numeric(fit0$phi) && length(fit0$phi) == 1L && fit0$phi > 0
)

expect("sigma2e is positive",
       is.numeric(fit0$sigma2e) && length(fit0$sigma2e) == 1L && fit0$sigma2e > 0
)

expect("sigma2b equals phi * sigma2e",
       abs(fit0$sigma2b - fit0$phi * fit0$sigma2e) < 1e-12
)

expect("method is exact_reml",
       identical(fit0$method, "exact_reml")
)

# =============================================================================
# 2. Fit with fixed effects
# =============================================================================
section("tune_reml: with X_fixed")

fit_xf <- tune_reml(
  y, A, Q_fun_simple,
  X_fixed = X,
  log_phi_lower = log(0.1),
  log_phi_upper = log(100),
  verbose = FALSE
)

expect("returns fastblm_tuned object with X_fixed",
       inherits(fit_xf, "fastblm_tuned")
)

expect("phi is positive with X_fixed",
       is.numeric(fit_xf$phi) && fit_xf$phi > 0
)

expect("sigma2e is positive with X_fixed",
       is.numeric(fit_xf$sigma2e) && fit_xf$sigma2e > 0
)

expect("n_eff equals n - rank(X)",
       identical(fit_xf$n_eff, n - qr(X)$rank)
)

# =============================================================================
# 3. Input validation
# =============================================================================
section("tune_reml: input validation")

expect_error("rejects X_fixed with wrong number of rows",
             tune_reml(
               y, A, Q_fun_simple,
               X_fixed = X[1:5, , drop = FALSE],
               verbose = FALSE
             )
)

expect_error("rejects unnamed theta_init",
             tune_reml(
               y, A, Q_fun_simple,
               theta_init = c(0),
               lower = -1,
               upper = 1,
               verbose = FALSE
             )
)

expect_error("rejects invalid logdet_method",
             tune_reml(
               y, A, Q_fun_simple,
               logdet_method = "banana",
               verbose = FALSE
             )
)

# =============================================================================
# 4. Rank-deficient X_fixed should error
# =============================================================================
section("tune_reml: rank-deficient X_fixed")

X_rankdef <- cbind(X, X[, 1])

expect_error("rank-deficient X_fixed errors",
             tune_reml(
               y, A, Q_fun_simple,
               X_fixed = X_rankdef,
               log_phi_lower = log(0.1),
               log_phi_upper = log(100),
               verbose = FALSE
             )
)

# =============================================================================
# 5. One-theta optimisation
# =============================================================================
section("tune_reml: one theta parameter")

Q_fun_theta <- function(theta) {
  scale <- exp(theta[["log_scale"]])
  list(
    Q = scale * Q_mat,
    Q_matrix = scale * Q_mat,
    log_det_Q = p * log(scale)
  )
}

fit_theta <- tune_reml(
  y, A, Q_fun_theta,
  theta_init = c(log_scale = 0),
  lower = -3,
  upper = 3,
  log_phi_lower = log(0.1),
  log_phi_upper = log(100),
  verbose = FALSE
)

expect("theta optimisation returns named theta",
       !is.null(names(fit_theta$theta)) && "log_scale" %in% names(fit_theta$theta)
)

expect("optimised theta is within bounds",
       fit_theta$theta[["log_scale"]] >= -3 &&
         fit_theta$theta[["log_scale"]] <= 3
)

expect("history is non-empty after theta optimisation",
       is.data.frame(fit_theta$history) && nrow(fit_theta$history) > 0
)

# =============================================================================
# 6. Constraint + fixed effects
# =============================================================================
section("tune_reml: constraint_matrix + X_fixed")

C_mat <- matrix(1, nrow = 1, ncol = p)

fit_constrained <- tune_reml(
  y, A, Q_fun_simple,
  X_fixed = X,
  constraint_matrix = C_mat,
  log_phi_lower = log(0.1),
  log_phi_upper = log(100),
  verbose = FALSE
)

expect("constrained + X_fixed fit returns fastblm_tuned",
       inherits(fit_constrained, "fastblm_tuned")
)

expect("constrained + X_fixed phi is positive",
       fit_constrained$phi > 0
)

expect("constrained + X_fixed sigma2e is positive",
       fit_constrained$sigma2e > 0
)

# =============================================================================
# 7. Exact objective matches dense reference at fixed phi
# =============================================================================
section("fixed-phi objective: dense reference agreement")

set.seed(123)
n_d <- 6
p_d <- 4

A_d <- matrix(rnorm(n_d * p_d), n_d, p_d)

B_d <- matrix(rnorm(p_d * p_d), p_d, p_d)
Q_d <- crossprod(B_d) + diag(0.5, p_d)

X_d <- cbind(1, seq_len(n_d))
R_d <- diag(c(1.0, 1.2, 0.8, 1.1, 0.9, 1.3))
Rinv_d <- solve(R_d)
C_d <- matrix(c(1, 1, 0, 0), nrow = 1)

phi_d <- 2.0
sigma2e_d <- 1.5

x_d <- drop(chol(sigma2e_d * phi_d * solve(Q_d)) %*% rnorm(p_d))
beta_d <- c(0.3, -0.2)
e_d <- drop(chol(sigma2e_d * R_d) %*% rnorm(n_d))
y_d <- drop(A_d %*% x_d + X_d %*% beta_d + e_d)

Q_fun_dense <- function(theta) {
  list(
    Q = Q_d,
    Q_matrix = Q_d,
    log_det_Q = as.numeric(determinant(Q_d, logarithm = TRUE)$modulus)
  )
}

cases <- list(
  list(name = "no_fixed_no_constraint", X = NULL, C = NULL),
  list(name = "fixed_no_constraint",    X = X_d,  C = NULL),
  list(name = "no_fixed_constraint",    X = NULL, C = C_d),
  list(name = "fixed_constraint",       X = X_d,  C = C_d)
)

dense_comp <- lapply(cases, function(case) {
  dense <- dense_reml_ll(
    y = y_d,
    A = A_d,
    Q = Q_d,
    phi = phi_d,
    X_fixed = case$X,
    R = R_d,
    C = case$C
  )

  new <- eval_new_ll(
    y = y_d,
    A = A_d,
    Q_fun = Q_fun_dense,
    phi = phi_d,
    X_fixed = case$X,
    R_inv = Rinv_d,
    C = case$C
  )

  data.frame(
    case = case$name,
    abs_diff = abs(dense$ll_drop_R_const - new$ll),
    sigma2e_diff = abs(dense$quad / dense$n_eff - new$sigma2e)
  )
})

dense_comp_df <- do.call(rbind, dense_comp)
print(dense_comp_df)

expect("fixed-phi ll matches dense reference in all four cases",
       max(dense_comp_df$abs_diff) < 1e-5
)

expect("fixed-phi sigma2e matches dense reference in all four cases",
       max(dense_comp_df$sigma2e_diff) < 1e-5
)

# =============================================================================
# 8. Constraint path agrees with dense reparameterisation
# =============================================================================
section("constraint_matrix: Schur correction vs reparameterisation")

set.seed(99)
n_c  <- 25
p_c  <- 6
r_c  <- 2
p_eff_c <- p_c - r_c

A_c <- matrix(rnorm(n_c * p_c), n_c, p_c)

Q_c_mat <- diag(p_c) + 0.2 * crossprod(matrix(rnorm(p_c^2), p_c, p_c)) / p_c
Q_c_mat <- (Q_c_mat + t(Q_c_mat)) / 2 + p_c * diag(p_c)
Q_c_mat <- Matrix::Matrix(Q_c_mat, sparse = FALSE)

logdet_Q_c <- as.numeric(Matrix::determinant(Q_c_mat, logarithm = TRUE)$modulus)
Q_fun_c <- function(theta) {
  list(
    Q = Q_c_mat,
    Q_matrix = Q_c_mat,
    log_det_Q = logdet_Q_c
  )
}

C_c <- matrix(rnorm(r_c * p_c), r_c, p_c)
while (qr(C_c)$rank < r_c) {
  C_c <- matrix(rnorm(r_c * p_c), r_c, p_c)
}

qr_Cc <- qr(t(C_c))
N_c <- qr.Q(qr_Cc, complete = TRUE)[, (r_c + 1L):p_c, drop = FALSE]

x_c <- N_c %*% rnorm(p_eff_c)
y_c <- as.numeric(A_c %*% x_c + rnorm(n_c, sd = sqrt(0.5)))

fit_schur <- tune_reml(
  y_c, A_c, Q_fun_c,
  constraint_matrix = C_c,
  log_phi_lower = log(0.1),
  log_phi_upper = log(50),
  verbose = FALSE
)

A_tilde_c <- A_c %*% N_c
Q_red_c   <- t(N_c) %*% as.matrix(Q_c_mat) %*% N_c
logdet_Qr <- as.numeric(determinant(Q_red_c, logarithm = TRUE)$modulus)

Q_fun_r <- function(theta) {
  list(
    Q = Matrix::Matrix(Q_red_c, sparse = FALSE),
    Q_matrix = Matrix::Matrix(Q_red_c, sparse = FALSE),
    log_det_Q = logdet_Qr
  )
}

fit_reparam <- tune_reml(
  y_c, A_tilde_c, Q_fun_r,
  log_phi_lower = log(0.1),
  log_phi_upper = log(50),
  verbose = FALSE
)

expect("Schur correction returns fastblm_tuned",
       inherits(fit_schur, "fastblm_tuned")
)

expect("Schur phi agrees with reparameterisation within 1%",
       abs(fit_schur$phi - fit_reparam$phi) / fit_reparam$phi < 0.01
)

expect("Schur sigma2e agrees with reparameterisation within 1%",
       abs(fit_schur$sigma2e - fit_reparam$sigma2e) / fit_reparam$sigma2e < 0.01
)

expect("Schur ll agrees with reparameterisation within 0.01",
       abs(fit_schur$value - fit_reparam$value) < 0.01
)

expect("Schur fit is finite and converged",
       is.finite(fit_schur$value) && fit_schur$optim$convergence == 0L
)

# =============================================================================
# 9. Print method
# =============================================================================
section("print.fastblm_tuned")

expect("print method runs without error",
       {
         tmp <- capture.output(print(fit0))
         length(tmp) > 0
       }
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
