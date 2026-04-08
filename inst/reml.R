#' Project out fixed effects for REML
#'
#' Projects y and A onto the orthogonal complement of X_fixed.
#' Only correct when R = I. For general R, whiten first.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param X_fixed n x q matrix of fixed covariates
#' @return list with projected y, A, and n_eff
#' @export
reml_project <- function(y, A, X_fixed) {
  if (!is_matrix(X_fixed)) stop("`X_fixed` must be a matrix.")
  if (nrow(X_fixed) != length(y)) stop("nrow(X_fixed) must equal length(y).")

  qr_fixed <- qr(as.matrix(X_fixed))
  q_rank   <- qr_fixed$rank

  if (q_rank < ncol(X_fixed))
    warning(sprintf("X_fixed has rank %d < ncol = %d.", q_rank, ncol(X_fixed)))

  list(
    y     = qr.resid(qr_fixed, y),
    A     = qr.resid(qr_fixed, as.matrix(A)),
    n_eff = length(y) - q_rank
  )
}

#' Tune hyperparameters via REML marginal likelihood
#'
#' Maximizes the REML marginal likelihood over theta (shape parameters) and
#' profiles out phi (signal-to-noise) and sigma2e analytically at each theta.
#'
#' When constraint_matrix C is supplied, the constrained marginal likelihood
#' p(y | Cx = 0) is computed via the Schur complement identity. The correction
#' costs r extra PCG solves per likelihood evaluation where r = rank(C).
#'
#' The concentrated constrained likelihood is:
#'   ll = -n_eff/2 * log(sigma2e) - 1/2*log|K| - p/2*log(phi)
#'        + 1/2*log|Q| - 1/2*log|M|
#' where M = I_r - Vt * H^{-1} * U  (r x r),
#'   U  = phi * A * Q^{-1} * C' * (C*Q^{-1}*C')^{-1}  (n x r)
#'   Vt = C * Q^{-1} * A'                               (r x n)
#'   H^{-1} v = v - A * K_H^{-1} * A' * v,  K_H = A'A + (1/phi)Q
#' and sigma2e = y' * H_c^{-1} * y / n_eff (profiled from constrained model).
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param Q_fun function(theta) -> list with Q, log_det_Q, Q_matrix
#' @param X_fixed n x q matrix of fixed-effect covariates (optional)
#' @param R_inv n x n inverse noise covariance. NULL = identity.
#' @param theta_init named numeric vector of initial theta values
#' @param lower lower bounds on theta
#' @param upper upper bounds on theta
#' @param log_phi_lower lower bound for phi search (log scale)
#' @param log_phi_upper upper bound for phi search (log scale)
#' @param logdet_method "lanczos" (default) or "cholesky"
#' @param n_lanczos_probes number of Lanczos probes
#' @param n_lanczos_steps number of Lanczos steps per probe
#' @param pcg_tol PCG tolerance
#' @param pcg_maxit max PCG iterations
#' @param constraint_matrix r x p matrix C such that C*x = 0
#' @param verbose logical
#' @return object of class fastblm_tuned
#' @export
tune_reml <- function(y, A, Q_fun,
                      X_fixed          = NULL,
                      R_inv            = NULL,
                      theta_init       = numeric(0),
                      lower            = rep(-Inf, length(theta_init)),
                      upper            = rep( Inf, length(theta_init)),
                      log_phi_lower    = log(0.01),
                      log_phi_upper    = log(1000),
                      logdet_method    = "lanczos",
                      n_lanczos_probes = 50L,
                      n_lanczos_steps  = 50L,
                      pcg_tol          = 1e-6,
                      pcg_maxit        = NULL,
                      constraint_matrix = NULL,
                      verbose          = TRUE) {

  y <- as.numeric(y)

  # --- project out fixed effects (REML) ------------------------------------
  if (!is.null(X_fixed)) {
    proj  <- reml_project(y, A, X_fixed)
    y     <- proj$y
    A     <- proj$A
    n_eff <- proj$n_eff
    if (verbose)
      message(sprintf("Projected out %d fixed effect(s); n_eff = %d.",
                      ncol(as.matrix(X_fixed)), n_eff))
  } else {
    n_eff <- length(y)
  }

  n <- length(y)
  p <- ncol(A)
  if (is.null(pcg_maxit)) pcg_maxit <- 4L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")

  # --- precompute once ------------------------------------------------------
  Rinv    <- resolve_Rinv(R_inv, n)
  Rinvy   <- as.numeric(Rinv %*% y)
  AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
  yRinvy  <- as.numeric(crossprod(y, Rinvy))

  AtRinvA       <- Matrix::crossprod(A, Rinv %*% A)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)

  .enrich_prior <- function(prior) {
    prior$AtRinvA_matrix <- AtRinvA
    prior
  }

  # --- validate constraint --------------------------------------------------
  # Build a closure that captures A, C_mat, pcg_tol, pcg_maxit correctly.
  # The closure is called inside .eval_reml_ll at each (phi, theta) evaluation.
  # It returns list(sigma2e_c, logdet_M) using the Schur complement identity.
  C_mat            <- NULL
  .eval_constraint <- NULL
  if (!is.null(constraint_matrix)) {
    if (!is_matrix(constraint_matrix))
      stop("`constraint_matrix` must be a matrix.")
    if (ncol(constraint_matrix) != p)
      stop(sprintf("`constraint_matrix` must have ncol = p = %d.", p))
    C_mat  <- as.matrix(constraint_matrix)
    r_constr <- qr(C_mat)$rank
    if (verbose)
      message(sprintf("Constraint: %d linear restriction(s) on x.", r_constr))

    # Capture A, C_mat, r in closure -- these never change across iterations
    A_dense_constr <- as.matrix(A)
    C_mat_constr   <- C_mat
    r_cap          <- r_constr
    pcg_tol_cap    <- pcg_tol
    pcg_maxit_cap  <- pcg_maxit

    y_constr <- y   # capture y for use inside closure
    .eval_constraint <- function(phi, apply_Q, prior, n_eff) {
      y <- y_constr
      # Requires Q as explicit matrix (for Q^{-1}C' solve)
      Q_mat <- if (!is.null(prior$Q_matrix)) prior$Q_matrix else
        if (is_matrix(prior$Q)) prior$Q else NULL
      if (is.null(Q_mat))
        stop("constraint_matrix requires Q to be available as an explicit matrix.")

      Q_mat_dense <- as.matrix(Q_mat)
      Q_inv_Ct    <- solve(Q_mat_dense, t(C_mat_constr))        # p x r
      CQinvCt     <- C_mat_constr %*% Q_inv_Ct                  # r x r
      CQinvCt_inv <- solve(CQinvCt)                             # r x r

      # U (n x r) and Vt (r x n)
      U  <- phi * A_dense_constr %*% (Q_inv_Ct %*% CQinvCt_inv)
      Vt <- t(Q_inv_Ct) %*% t(A_dense_constr)

      # K_H = A'A + (1/phi)Q  (no sigma2e); H^{-1} = I - A*K_H^{-1}*A'
      apply_AtRinvA_loc <- function(v) as.numeric(
        crossprod(A_dense_constr, A_dense_constr %*% v))  # A'Av (R=I)
      apply_K_H <- function(v) apply_AtRinvA_loc(v) + (1/phi) * apply_Q(v)

      Hinv_apply <- function(V) {
        K_H_inv_AtV <- apply(
          crossprod(A_dense_constr, V), 2,
          function(b) pcg(apply_K_H, b, tol = pcg_tol_cap,
                          maxit = pcg_maxit_cap)$x
        )
        V - A_dense_constr %*% K_H_inv_AtV
      }

      HinvU <- Hinv_apply(U)                                    # n x r
      M     <- diag(r_cap) - Vt %*% HinvU                      # r x r
      Minv  <- solve(M)

      Hinvy <- Hinv_apply(matrix(y, ncol = 1L))[, 1L]
      s_c   <- as.numeric(crossprod(y,
                                    Hinvy + HinvU %*% Minv %*% (Vt %*% Hinvy)))
      if (!is.finite(s_c) || s_c <= 0)
        return(list(sigma2e = NA, logdet_M = NA))

      logdet_M <- as.numeric(determinant(M, logarithm = TRUE)$modulus)
      list(sigma2e = s_c / n_eff, logdet_M = logdet_M)
    }
  }

  # pre-draw Lanczos probes once for deterministic objective
  set.seed(1L)
  probes <- matrix(sample(c(-1L, 1L), p * n_lanczos_probes, replace = TRUE),
                   nrow = p, ncol = n_lanczos_probes)

  x_warm <- rep(0, p)

  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta -------------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)
    prior <- .enrich_prior(prior)

    apply_Q <- as_apply(prior$Q)

    logdet_Q <- if (!is.null(prior$log_det_Q)) {
      prior$log_det_Q
    } else {
      lanczos_logdet(apply_Q, probes, n_lanczos_steps)
    }
    if (!is.finite(logdet_Q)) return(.Machine$double.xmax)

    phi_hat <- .profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
      log_phi_lower, log_phi_upper,
      logdet_method, pcg_tol, pcg_maxit, x_warm,
      eval_constraint = .eval_constraint
    )

    ll_res <- .eval_reml_ll(
      phi_hat, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
      logdet_method, pcg_tol, pcg_maxit, x_warm,
      eval_constraint = .eval_constraint
    )
    if (!is.finite(ll_res$ll)) return(.Machine$double.xmax)

    x_warm <<- ll_res$x

    if (verbose) {
      theta_str <- if (length(theta) > 0L)
        paste(sprintf("%s=%.4g", names(theta), theta), collapse = "  ")
      else ""
      message(sprintf("  iter %d: %s  log_phi=%.3f  sigma2e=%.4g  ll=%.4f",
                      iter_count, theta_str, log(phi_hat),
                      ll_res$sigma2e, ll_res$ll))
    }

    history[[iter_count]] <<- list(
      theta   = theta,
      phi     = phi_hat,
      sigma2e = ll_res$sigma2e,
      ll      = ll_res$ll
    )

    -ll_res$ll
  }

  # --- dispatch optimizer ---------------------------------------------------
  if (verbose) message("Starting REML optimisation...")

  if (length(theta_init) == 0L) {
    prior_0   <- .enrich_prior(Q_fun(numeric(0)))
    apply_Q0  <- as_apply(prior_0$Q)
    logdet_Q0 <- prior_0$log_det_Q %||%
      lanczos_logdet(apply_Q0, probes, n_lanczos_steps)

    phi_opt <- stats::optimize(
      f        = function(lp) {
        res <- .eval_reml_ll(
          exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q0, prior_0,
          p, n, n_eff, logdet_Q0, probes, n_lanczos_steps,
          logdet_method, pcg_tol, pcg_maxit, x_warm,
          eval_constraint = .eval_constraint
        )
        -res$ll
      },
      interval = c(log_phi_lower, log_phi_upper)
    )
    theta_opt <- numeric(0)
    phi_opt   <- exp(phi_opt$minimum)
    optim_res <- list(convergence = 0L)

  } else if (length(theta_init) == 1L) {
    opt <- stats::optimize(
      f        = function(th) {
        theta <- th; names(theta) <- names(theta_init)
        .objective(theta)
      },
      interval = c(unname(lower), unname(upper))
    )
    theta_opt        <- opt$minimum
    names(theta_opt) <- names(theta_init)
    optim_res        <- list(convergence = 0L)

  } else {
    opt <- stats::optim(
      par     = unname(theta_init),
      fn      = function(par) {
        theta <- par; names(theta) <- names(theta_init)
        .objective(theta)
      },
      method  = "L-BFGS-B",
      lower   = unname(lower),
      upper   = unname(upper),
      control = list(maxit = 300L, trace = if (verbose) 1L else 0L)
    )
    theta_opt        <- opt$par
    names(theta_opt) <- names(theta_init)
    optim_res        <- opt
  }

  # --- recover at optimum ---------------------------------------------------
  prior_opt    <- .enrich_prior(Q_fun(theta_opt))
  apply_Q_opt  <- as_apply(prior_opt$Q)
  logdet_Q_opt <- prior_opt$log_det_Q %||%
    lanczos_logdet(apply_Q_opt, probes, n_lanczos_steps)

  phi_final <- .profile_phi(
    AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, n_eff, logdet_Q_opt, probes, n_lanczos_steps,
    log_phi_lower, log_phi_upper,
    logdet_method, pcg_tol, pcg_maxit, x_warm,
    eval_constraint = .eval_constraint
  )

  ll_final <- .eval_reml_ll(
    phi_final, AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, n_eff, logdet_Q_opt, probes, n_lanczos_steps,
    logdet_method, pcg_tol, pcg_maxit, x_warm,
    eval_constraint = .eval_constraint
  )

  if (verbose)
    message(sprintf("Optimum: phi=%.4g  sigma2e=%.4g  ll=%.4f",
                    phi_final, ll_final$sigma2e, ll_final$ll))

  hist_df <- if (length(history) > 0L)
    do.call(rbind, lapply(history, as.data.frame))
  else data.frame()

  structure(
    list(
      theta   = theta_opt,
      phi     = phi_final,
      sigma2e = ll_final$sigma2e,
      sigma2b = phi_final * ll_final$sigma2e,
      Q       = prior_opt$Q,
      value   = ll_final$ll,
      method  = "reml",
      optim   = optim_res,
      history = hist_df
    ),
    class = "fastblm_tuned"
  )
}

# -----------------------------------------------------------------------
# Internal: profile phi via 1D golden section at fixed theta
# -----------------------------------------------------------------------
.profile_phi <- function(AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
                         p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
                         log_phi_lower, log_phi_upper,
                         logdet_method, pcg_tol, pcg_maxit, x_warm,
                         eval_constraint = NULL) {
  opt <- stats::optimize(
    f = function(lp) {
      res <- .eval_reml_ll(
        exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
        logdet_method, pcg_tol, pcg_maxit, x_warm,
        eval_constraint = eval_constraint
      )
      -res$ll
    },
    interval = c(log_phi_lower, log_phi_upper),
    tol      = 1e-4
  )
  exp(opt$minimum)
}

# -----------------------------------------------------------------------
# Internal: evaluate REML log likelihood at fixed (phi, theta)
#
# Unconstrained (C_mat = NULL):
#   ll = -n_eff/2 * log(sigma2e) - 1/2*log|K| - p/2*log(phi) + 1/2*log|Q|
#   sigma2e = (y'R^{-1}y - x'A'R^{-1}y) / n_eff
#
# Constrained (C_mat supplied):
#   ll = -n_eff/2 * log(sigma2e_c) - 1/2*log|K| - p/2*log(phi)
#        + 1/2*log|Q| - 1/2*log|M|
#   sigma2e_c = y' H_c^{-1} y / n_eff
#   M = I_r - Vt * H^{-1} * U  (r x r)
#   H^{-1} v = v - A * K_H^{-1} * A' * v  (K_H = A'A + (1/phi)Q, via PCG)
#
# K = A'R^{-1}A + (1/phi)*Q  (posterior precision)
# n_eff = n - rank(X_fixed)
# -----------------------------------------------------------------------
.eval_reml_ll <- function(phi, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
                          p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
                          logdet_method, pcg_tol, pcg_maxit, x_warm,
                          eval_constraint = NULL) {

  apply_K <- function(v) apply_AtRinvA(v) + (1/phi) * apply_Q(v)

  pcg_res <- tryCatch(
    pcg(apply_K, AtRinvy, x0 = x_warm, tol = pcg_tol, maxit = pcg_maxit),
    error = function(e) list(converged = FALSE)
  )
  if (!pcg_res$converged) return(list(ll = -Inf, sigma2e = NA, x = x_warm))
  x <- pcg_res$x

  # --- log|K| ---------------------------------------------------------------
  logdet_K <- .eval_logdet_K(
    apply_K, prior, phi, p, probes, n_lanczos_steps, logdet_method
  )
  if (!is.finite(logdet_K)) return(list(ll = -Inf, sigma2e = NA, x = x))

  if (is.null(eval_constraint)) {
    # --- Unconstrained -------------------------------------------------------
    yHinvy  <- yRinvy - as.numeric(crossprod(AtRinvy, x))
    if (yHinvy <= 0) return(list(ll = -Inf, sigma2e = NA, x = x))
    sigma2e <- yHinvy / n_eff

    ll <- -n_eff/2 * log(sigma2e) -
      1/2     * logdet_K     -
      p/2     * log(phi)     +
      1/2     * logdet_Q

  } else {
    # --- Constrained: Schur correction ---------------------------------------
    # Delegate to the closure built in tune_reml which captures A correctly.
    constr <- eval_constraint(phi, apply_Q, prior, n_eff)
    if (is.na(constr$sigma2e) || !is.finite(constr$logdet_M))
      return(list(ll = -Inf, sigma2e = NA, x = x))

    sigma2e  <- constr$sigma2e
    logdet_M <- constr$logdet_M

    ll <- -n_eff/2 * log(sigma2e) -
      1/2     * logdet_K     -
      p/2     * log(phi)     +
      1/2     * logdet_Q     -
      1/2     * logdet_M
  }

  list(ll = ll, sigma2e = sigma2e, x = x)
}

# -----------------------------------------------------------------------
# Internal: log|K| via Lanczos or Cholesky
# -----------------------------------------------------------------------
.eval_logdet_K <- function(apply_K, prior, phi, p,
                           probes, n_lanczos_steps, logdet_method) {

  if (logdet_method == "cholesky") {
    Q_mat <- if (!is.null(prior$Q_matrix)) {
      prior$Q_matrix
    } else if (is_matrix(prior$Q)) {
      prior$Q
    } else {
      message("logdet_method='cholesky' requires Q_matrix. Falling back to Lanczos.")
      return(lanczos_logdet(apply_K, probes, n_lanczos_steps))
    }
    if (!is.null(prior$AtRinvA_matrix)) {
      K_mat <- Matrix::forceSymmetric(prior$AtRinvA_matrix + (1/phi) * Q_mat)
      tryCatch(
        as.numeric(Matrix::determinant(K_mat, logarithm = TRUE)$modulus),
        error = function(e) lanczos_logdet(apply_K, probes, n_lanczos_steps)
      )
    } else {
      lanczos_logdet(apply_K, probes, n_lanczos_steps)
    }
  } else {
    lanczos_logdet(apply_K, probes, n_lanczos_steps)
  }
}

# -----------------------------------------------------------------------
# Print method
# -----------------------------------------------------------------------
#' @export
print.fastblm_tuned <- function(x, ...) {
  cat("fastblm_tuned\n")
  cat(sprintf("  method  : %s\n", x$method))
  if (length(x$theta) > 0L) {
    cat("  theta   :\n")
    for (nm in names(x$theta))
      cat(sprintf("    %s = %.6g\n", nm, x$theta[[nm]]))
  }
  cat(sprintf("  phi     : %.6g\n", x$phi))
  cat(sprintf("  sigma2e : %.6g\n", x$sigma2e))
  cat(sprintf("  sigma2b : %.6g\n", x$sigma2b))
  cat(sprintf("  ll      : %.4f\n", x$value))
  cat(sprintf("  converge: %d\n",   x$optim$convergence))
  invisible(x)
}
