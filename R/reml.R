#' Project out fixed effects for REML
#'
#' Projects y and A onto the orthogonal complement of X_fixed.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param X_fixed n x q matrix of fixed covariates
#'
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
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param Q_fun function(theta) -> list with:
#'   \itemize{
#'     \item Q: sparse matrix or function v -> Qv (required)
#'     \item log_det_Q: scalar log|Q| (optional, triggers Lanczos if missing)
#'     \item Q_matrix: explicit sparse matrix (optional, needed for exact logdet_K)
#'   }
#' @param R_inv n x n inverse noise covariance. NULL = identity.
#' @param theta_init named numeric vector of initial theta values
#' @param lower lower bounds on theta
#' @param upper upper bounds on theta
#' @param log_phi_lower lower bound for phi search (log scale)
#' @param log_phi_upper upper bound for phi search (log scale)
#' @param logdet_method "lanczos" (default) or "cholesky"
#' @param n_lanczos_probes number of Lanczos probes
#' @param n_lanczos_steps number of Lanczos steps per probe
#' @param pcg_tol PCG tolerance for linear solves
#' @param pcg_maxit max PCG iterations
#' @param verbose logical
#'
#' @return object of class fastblm_tuned
#' @export
tune_reml <- function(y, A, Q_fun,
                      R_inv           = NULL,
                      theta_init      = numeric(0),
                      lower           = rep(-Inf, length(theta_init)),
                      upper           = rep( Inf, length(theta_init)),
                      log_phi_lower   = log(0.01),
                      log_phi_upper   = log(1000),
                      logdet_method   = "lanczos",
                      n_lanczos_probes = 50L,
                      n_lanczos_steps  = 50L,
                      pcg_tol         = 1e-6,
                      pcg_maxit       = NULL,
                      verbose         = TRUE) {

  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(A)
  if (is.null(pcg_maxit)) pcg_maxit <- 4L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")

  # --- precompute once ---------------------------------------------------
  Rinv    <- resolve_Rinv(R_inv, n)
  Rinvy   <- as.numeric(Rinv %*% y)
  AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
  yRinvy  <- as.numeric(crossprod(y, Rinvy))

  # precompute A'R^{-1}A once -- this is O(n*p^2) and must not be repeated
  # store as sparse matrix if possible, otherwise as operator
  AtRinvA       <- Matrix::crossprod(A, Rinv %*% A)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)

  # helper to inject AtRinvA into prior for Cholesky logdet path
  .enrich_prior <- function(prior) {
    prior$AtRinvA_matrix <- AtRinvA
    prior
  }

  # pre-draw Lanczos probes once for deterministic objective
  set.seed(1L)
  probes <- matrix(sample(c(-1L, 1L), p * n_lanczos_probes, replace = TRUE),
                   nrow = p, ncol = n_lanczos_probes)

  # warm start for PCG across iterations
  x_warm <- rep(0, p)

  # --- iteration history -------------------------------------------------
  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta ----------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)
    prior <- .enrich_prior(prior)

    apply_Q  <- as_apply(prior$Q)

    # log|Q| -- use supplied value or Lanczos
    logdet_Q <- if (!is.null(prior$log_det_Q)) {
      prior$log_det_Q
    } else {
      lanczos_logdet(apply_Q, probes, n_lanczos_steps)
    }
    if (!is.finite(logdet_Q)) return(.Machine$double.xmax)

    # profile phi at this theta
    phi_hat <- .profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdet_Q, probes, n_lanczos_steps,
      log_phi_lower, log_phi_upper,
      logdet_method, pcg_tol, pcg_maxit, x_warm
    )

    # evaluate ll at phi_hat
    ll_res <- .eval_reml_ll(
      phi_hat, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, logdet_Q, probes, n_lanczos_steps,
      logdet_method, pcg_tol, pcg_maxit, x_warm
    )
    if (!is.finite(ll_res$ll)) return(.Machine$double.xmax)

    # update warm start
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

  # --- dispatch optimizer ------------------------------------------------
  if (verbose) message("Starting REML optimisation...")

  if (length(theta_init) == 0L) {
    # no shape parameters -- just profile phi
    prior_0  <- .enrich_prior(Q_fun(numeric(0)))
    apply_Q0 <- as_apply(prior_0$Q)
    logdet_Q0 <- prior_0$log_det_Q %||%
      lanczos_logdet(apply_Q0, probes, n_lanczos_steps)

    phi_opt <- stats::optimize(
      f        = function(lp) {
        res <- .eval_reml_ll(
          exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q0, prior_0,
          p, n, logdet_Q0, probes, n_lanczos_steps,
          logdet_method, pcg_tol, pcg_maxit, x_warm
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

  # --- recover at optimum -----------------------------------------------
  prior_opt    <- .enrich_prior(Q_fun(theta_opt))
  apply_Q_opt  <- as_apply(prior_opt$Q)
  logdet_Q_opt <- prior_opt$log_det_Q %||%
    lanczos_logdet(apply_Q_opt, probes, n_lanczos_steps)

  phi_final <- .profile_phi(
    AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, logdet_Q_opt, probes, n_lanczos_steps,
    log_phi_lower, log_phi_upper,
    logdet_method, pcg_tol, pcg_maxit, x_warm
  )

  ll_final <- .eval_reml_ll(
    phi_final, AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, logdet_Q_opt, probes, n_lanczos_steps,
    logdet_method, pcg_tol, pcg_maxit, x_warm
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
                         p, n, logdet_Q, probes, n_lanczos_steps,
                         log_phi_lower, log_phi_upper,
                         logdet_method, pcg_tol, pcg_maxit, x_warm) {
  opt <- stats::optimize(
    f = function(lp) {
      res <- .eval_reml_ll(
        exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n, logdet_Q, probes, n_lanczos_steps,
        logdet_method, pcg_tol, pcg_maxit, x_warm
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
# ll = -n/2 * log(sigma2e) - 1/2 * log|K| - p/2 * log(phi) + 1/2 * log|Q|
#
# where K = A'R^{-1}A + (1/phi) Q
# and sigma2e = (y'R^{-1}y - x' A'R^{-1}y) / n  (profiled out)
# -----------------------------------------------------------------------
.eval_reml_ll <- function(phi, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
                          p, n, logdet_Q, probes, n_lanczos_steps,
                          logdet_method, pcg_tol, pcg_maxit, x_warm) {

  # apply_K = A'R^{-1}A v + (1/phi) Q v
  # apply_AtRinvA is precomputed once outside the loop
  apply_K <- function(v) apply_AtRinvA(v) + (1/phi) * apply_Q(v)

  pcg_res <- tryCatch(
    pcg(apply_K, AtRinvy, x0 = x_warm, tol = pcg_tol, maxit = pcg_maxit),
    error = function(e) list(converged = FALSE)
  )
  if (!pcg_res$converged) return(list(ll = -Inf, sigma2e = NA, x = x_warm))
  x <- pcg_res$x

  # --- profile sigma2e --------------------------------------------------
  yHinvy  <- yRinvy - as.numeric(crossprod(AtRinvy, x))
  if (yHinvy <= 0) return(list(ll = -Inf, sigma2e = NA, x = x))
  sigma2e <- yHinvy / n

  # --- log|K| -----------------------------------------------------------
  logdet_K <- .eval_logdet_K(
    apply_K, prior, phi, p,
    probes, n_lanczos_steps, logdet_method
  )
  if (!is.finite(logdet_K)) return(list(ll = -Inf, sigma2e = NA, x = x))

  # --- REML ll ----------------------------------------------------------
  ll <- -n/2 * log(sigma2e) -
    1/2  * logdet_K     -
    p/2  * log(phi)     +
    1/2  * logdet_Q

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
      message("logdet_method = 'cholesky' requires Q_matrix in prior. Falling back to Lanczos.")
      return(lanczos_logdet(apply_K, probes, n_lanczos_steps))
    }
    # need AtRinvA as matrix too -- stored in prior if supplied
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
