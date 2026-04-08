#' Tune hyperparameters via REML marginal likelihood
#'
#' Maximizes the REML marginal likelihood over theta (shape parameters) and
#' profiles out phi (signal-to-noise) and sigma2e analytically at each theta.
#'
#' Fixed effects are handled by implicit projection -- the annihilator
#' \eqn{M_X = I - X(X^\top X)^{-1}X^\top} is applied to y and to A matvecs
#' without ever forming the dense projected matrix M_X A. This keeps A sparse
#' throughout.
#'
#' Three logdet methods are available:
#' \describe{
#'   \item{\code{"woodbury"}}{Exact. Uses the matrix determinant lemma:
#'     \eqn{\log|K_{\text{proj}}| = \log|(1/\phi)Q| + \log|I + \phi B|}
#'     where \eqn{B = M_X A Q^{-1} A^\top M_X^\top} is n x n. B is formed
#'     once per theta evaluation (n Q^{-1} solves), then \code{chol(I+phi*B)}
#'     is computed per phi evaluation (~1s for n=7000). Exact and correct
#'     when X_fixed is present -- the only method that accounts for the
#'     near-collinearity of X with the column space of A.}
#'   \item{\code{"cholesky"}}{Exact. Forms K = A'R^{-1}A + (1/phi)Q explicitly
#'     and takes its sparse Cholesky. Fast when A'A is sparse and p is moderate.
#'     Does NOT correctly account for fixed-effect projection -- use
#'     \code{"woodbury"} when X_fixed is supplied.}
#'   \item{\code{"lanczos"}}{Stochastic. WARNING: severely biased when the
#'     matrix has a large near-null space (rank(A) << p), which is the typical
#'     downscaling regime. Use only when A'A is full rank or as a fast
#'     approximation where bias is acceptable.}
#' }
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix (sparse dgCMatrix recommended)
#' @param Q_fun function(theta) returning a list with elements:
#'   \code{Q} (p x p precision matrix or function),
#'   \code{Q_matrix} (explicit p x p sparse matrix, required for
#'   \code{logdet_method = "cholesky"}),
#'   \code{log_det_Q} (optional precomputed scalar).
#' @param X_fixed n x q matrix of fixed-effect covariates (optional).
#'   Fixed effects are profiled out via implicit M_X projection.
#'   Use \code{logdet_method = "woodbury"} when X_fixed is supplied.
#' @param R_inv n x n inverse noise covariance. NULL = identity.
#' @param theta_init named numeric vector of initial theta values
#' @param lower lower bounds on theta
#' @param upper upper bounds on theta
#' @param log_phi_lower lower bound for phi search (log scale)
#' @param log_phi_upper upper bound for phi search (log scale)
#' @param solver one of \code{"cholesky"}, \code{"pcg"}. Use \code{"pcg"}
#'   for large p. When \code{logdet_method = "woodbury"}, \code{"pcg"} is
#'   always used for the linear solve regardless of this argument.
#' @param logdet_method one of \code{"woodbury"}, \code{"cholesky"},
#'   \code{"lanczos"}. Default \code{"lanczos"}. See Details.
#' @param n_lanczos_probes Lanczos probes (only used when
#'   \code{logdet_method = "lanczos"})
#' @param n_lanczos_steps Lanczos steps per probe
#' @param pcg_tol PCG tolerance for tuning evaluations. Default 1e-3.
#' @param pcg_tol_final PCG tolerance for the final recover-at-optimum solve.
#' @param pcg_maxit max PCG iterations
#' @param precond optional PCG preconditioner function \code{v -> M^{-1}v}.
#' @param constraint_matrix r x p matrix C such that \eqn{Cx = 0}
#' @param verbose logical
#'
#' @return object of class fastblm_tuned
#' @export
tune_reml <- function(y, A, Q_fun,
                      X_fixed           = NULL,
                      R_inv             = NULL,
                      theta_init        = numeric(0),
                      lower             = rep(-Inf, length(theta_init)),
                      upper             = rep( Inf, length(theta_init)),
                      log_phi_lower     = log(0.01),
                      log_phi_upper     = log(1000),
                      solver            = c("cholesky", "pcg"),
                      logdet_method     = "lanczos",
                      n_lanczos_probes  = 50L,
                      n_lanczos_steps   = 50L,
                      pcg_tol           = 1e-3,
                      pcg_tol_final     = 1e-6,
                      pcg_maxit         = NULL,
                      precond           = NULL,
                      constraint_matrix = NULL,
                      verbose           = TRUE) {

  y      <- as.numeric(y)
  solver <- match.arg(solver)
  precond_apply <- if (!is.null(precond)) as_apply(precond) else NULL

  n <- length(y)
  p <- ncol(A)

  # Woodbury logdet forces PCG for the linear solve -- Cholesky of K is not
  # needed and would be inconsistent with the projected system
  if (logdet_method == "woodbury") solver <- "pcg"

  # -------------------------------------------------------------------------
  # Fixed-effect projection setup
  #
  # M_X = I - X(X'X)^{-1}X'  applied implicitly:
  #   y_proj   = M_X y               (n-vector, once)
  #   AtRinvy  = A'R^{-1}(M_X y)     (p-vector, once)
  #   operator: v -> A'R^{-1}(M_X Av)  (per PCG iteration, sparse)
  #
  # For logdet_method = "woodbury", B = M_X A Q^{-1} A' M_X' is formed
  # once per theta by applying Q^{-1} to each row of A then projecting.
  # logdet(K_proj) = logdet((1/phi)Q) + logdet(I + phi*B)  -- exact.
  # -------------------------------------------------------------------------

  if (!is.null(X_fixed)) {
    X_fixed <- as.matrix(X_fixed)
    q       <- ncol(X_fixed)
    n_eff   <- n - q
    XtX_inv <- solve(crossprod(X_fixed))
    MX      <- function(u) {
      if (is.matrix(u))
        u - X_fixed %*% (XtX_inv %*% crossprod(X_fixed, u))
      else
        u - as.numeric(X_fixed %*% (XtX_inv %*% crossprod(X_fixed, u)))
    }
    if (verbose)
      message(sprintf(
        "REML: projecting out %d fixed effect(s) implicitly, n_eff = %d.",
        q, n_eff))
  } else {
    n_eff <- n
    MX    <- NULL
    q     <- 0L
  }

  if (is.null(pcg_maxit)) pcg_maxit <- 4L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")

  # --- precompute fixed quantities ------------------------------------------
  Rinv  <- resolve_Rinv(R_inv, n)
  Rinvy <- as.numeric(Rinv %*% y)

  if (!is.null(MX)) {
    y_proj     <- MX(y)
    Rinvy_proj <- as.numeric(Rinv %*% y_proj)
    AtRinvy    <- as.numeric(Matrix::crossprod(A, Rinvy_proj))
    yRinvy     <- as.numeric(crossprod(y_proj, Rinvy_proj))
    AtRinvX    <- as.matrix(Matrix::crossprod(A, Rinv %*% X_fixed))  # p x q

    # Projected operator: v -> A'R^{-1} M_X A v  (never densifies A)
    apply_AtRinvA <- function(v) {
      Av      <- as.numeric(A %*% v)
      XtAv    <- as.numeric(crossprod(X_fixed, Av))
      RinvMAv <- as.numeric(Rinv %*% (Av - X_fixed %*% (XtX_inv %*% XtAv)))
      as.numeric(Matrix::crossprod(A, RinvMAv))
    }

    # For non-woodbury Cholesky logdet: need unprojected A'R^{-1}A
    AtRinvA <- if (logdet_method == "cholesky") {
      if (verbose) message("Forming A'R^{-1}A for Cholesky logdet...")
      Matrix::crossprod(A, Rinv %*% A)
    } else {
      NULL
    }

    # REML log-det correction: -1/2 log|X'H^{-1}X|
    # Only needed for non-woodbury paths; woodbury handles this exactly
    XtRinvX <- as.matrix(crossprod(X_fixed, Rinv %*% X_fixed))
    .logdet_correction <- if (logdet_method != "woodbury") {
      function(apply_K, p) {
        KinvAtRinvX <- vapply(seq_len(q), function(j)
          pcg(apply_K, AtRinvX[, j], tol = pcg_tol_final,
              maxit = pcg_maxit)$x,
          numeric(p))
        XtHinvX <- XtRinvX - crossprod(AtRinvX, KinvAtRinvX)
        ld <- tryCatch(
          as.numeric(determinant(XtHinvX, logarithm = TRUE)$modulus),
          error = function(e) NA_real_
        )
        if (!is.finite(ld)) return(NA_real_)
        -0.5 * ld
      }
    } else {
      NULL  # woodbury logdet is already exact and accounts for projection
    }

  } else {
    # No fixed effects
    AtRinvy <- as.numeric(Matrix::crossprod(A, Rinvy))
    yRinvy  <- as.numeric(crossprod(y, Rinvy))
    AtRinvA <- if (logdet_method == "cholesky") {
      if (verbose) message("Forming A'R^{-1}A for Cholesky logdet...")
      Matrix::crossprod(A, Rinv %*% A)
    } else {
      NULL
    }
    apply_AtRinvA      <- function(v)
      as.numeric(Matrix::crossprod(A, Rinv %*% (A %*% v)))
    AtRinvX            <- NULL
    .logdet_correction <- NULL
  }

  .enrich_prior <- function(prior) {
    prior$AtRinvA_matrix <- AtRinvA
    prior
  }

  # -------------------------------------------------------------------------
  # Woodbury logdet: form B = M_X A Q^{-1} A' M_X' once per theta
  #
  # logdet(K_proj) = -p*log(phi) + logdet(Q) + logdet(I + phi*B)
  #
  # B is n x n, formed via n Q^{-1} solves (~31s at n=6850, p=20000).
  # Then per-phi: chol(I + phi*B) costs ~1s for n=6850.
  #
  # When X_fixed is NULL, M_X = I and B = A Q^{-1} A'.
  # -------------------------------------------------------------------------
  .make_B <- function(CQ) {
    if (verbose) message("  Forming B = MX A Q^{-1} A' MX' (n Q^{-1} solves)...")
    t0 <- proc.time()["elapsed"]
    # Q^{-1} applied to each row of A (= each col of A')
    QinvAt <- vapply(seq_len(n), function(j) {
      aj <- as.numeric(A[j, , drop = FALSE])
      as.numeric(Matrix::solve(CQ, aj))
    }, numeric(p))                              # p x n
    AQinvAt <- as.matrix(A %*% QinvAt)         # n x n
    # Apply M_X on both sides if X_fixed present
    B <- if (!is.null(MX)) MX(t(MX(t(AQinvAt)))) else AQinvAt
    if (verbose)
      message(sprintf("  B formed in %.1fs", proc.time()["elapsed"] - t0))
    B
  }

  .logdet_K_woodbury <- function(phi, B, logdet_Q) {
    # logdet(K_proj) = -p*log(phi) + logdet(Q) + logdet(I + phi*B)
    CM_B <- tryCatch(chol(diag(n) + phi * B), error = function(e) NULL)
    if (is.null(CM_B)) return(-Inf)
    -p * log(phi) + logdet_Q + 2 * sum(log(diag(CM_B)))
  }

  # --- constraint -----------------------------------------------------------
  C_mat            <- NULL
  .eval_constraint <- NULL
  if (!is.null(constraint_matrix)) {
    if (!is_matrix(constraint_matrix))
      stop("`constraint_matrix` must be a matrix.")
    if (ncol(constraint_matrix) != p)
      stop(sprintf("`constraint_matrix` must have ncol = p = %d.", p))
    C_mat    <- as.matrix(constraint_matrix)
    r_constr <- qr(C_mat)$rank
    if (verbose)
      message(sprintf("Constraint: %d linear restriction(s) on x.", r_constr))

    A_dense_constr <- as.matrix(A)
    C_mat_constr   <- C_mat
    Rinv_constr    <- Rinv
    r_cap          <- r_constr
    pcg_tol_cap    <- pcg_tol
    pcg_maxit_cap  <- pcg_maxit
    y_constr       <- y

    .eval_constraint <- function(phi, apply_Q, prior, n_eff) {
      y <- y_constr
      Q_mat <- if (!is.null(prior$Q_matrix)) prior$Q_matrix else
        if (is_matrix(prior$Q)) prior$Q else NULL
      if (is.null(Q_mat))
        stop("constraint_matrix requires Q to be available as an explicit matrix.")

      Q_mat_dense <- as.matrix(Q_mat)
      Q_inv_Ct    <- solve(Q_mat_dense, t(C_mat_constr))
      CQinvCt     <- C_mat_constr %*% Q_inv_Ct
      CQinvCt_inv <- solve(CQinvCt)

      RinvA <- as.matrix(Rinv_constr %*% A_dense_constr)
      apply_AtRinvA_loc <- function(v)
        as.numeric(crossprod(A_dense_constr, RinvA %*% v))
      apply_K_loc <- function(v)
        apply_AtRinvA_loc(v) + (1/phi) * apply_Q(v)

      Hinv_apply <- function(V) {
        RinvV      <- as.matrix(Rinv_constr %*% V)
        AtRinvV    <- crossprod(A_dense_constr, RinvV)
        K_inv_AtRV <- apply(AtRinvV, 2,
                            function(b) pcg(apply_K_loc, b,
                                            tol   = pcg_tol_cap,
                                            maxit = pcg_maxit_cap)$x)
        RinvV - RinvA %*% K_inv_AtRV
      }

      U  <- phi * A_dense_constr %*% (Q_inv_Ct %*% CQinvCt_inv)
      Vt <- t(Q_inv_Ct) %*% t(RinvA)

      HinvU <- Hinv_apply(U)
      M     <- diag(r_cap) - Vt %*% HinvU
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

  set.seed(1L)
  probes <- matrix(sample(c(-1L, 1L), p * n_lanczos_probes, replace = TRUE),
                   nrow = p, ncol = n_lanczos_probes)

  x_warm     <- rep(0, p)
  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta -------------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)
    prior <- .enrich_prior(prior)

    apply_Q  <- as_apply(prior$Q)

    # logdet_Q: from prior if provided, else Lanczos (Q alone is well-behaved)
    logdet_Q <- if (!is.null(prior$log_det_Q)) prior$log_det_Q else
      lanczos_logdet(apply_Q, probes, n_lanczos_steps)
    if (!is.finite(logdet_Q)) return(.Machine$double.xmax)

    # For woodbury: form B once per theta, cache for phi profile
    B_cache <- if (logdet_method == "woodbury") {
      Q_mat <- if (!is.null(prior$Q_matrix)) prior$Q_matrix else
        if (is_matrix(prior$Q)) prior$Q else NULL
      if (is.null(Q_mat)) {
        warning("logdet_method='woodbury' requires Q_matrix in Q_fun output.")
        return(.Machine$double.xmax)
      }
      CQ <- tryCatch(
        Matrix::Cholesky(Q_mat, LDL = FALSE, perm = TRUE),
        error = function(e) NULL
      )
      if (is.null(CQ)) return(.Machine$double.xmax)
      .make_B(CQ)
    } else {
      NULL
    }

    phi_hat <- .profile_phi(
      AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
      log_phi_lower, log_phi_upper,
      logdet_method, pcg_tol, pcg_maxit, x_warm,
      eval_constraint   = .eval_constraint,
      logdet_correction = .logdet_correction,
      solver            = solver,
      precond           = precond_apply,
      B_cache           = B_cache
    )

    ll_res <- .eval_reml_ll(
      phi_hat, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
      p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
      logdet_method, pcg_tol, pcg_maxit, x_warm,
      eval_constraint   = .eval_constraint,
      logdet_correction = .logdet_correction,
      solver            = solver,
      precond           = precond_apply,
      B_cache           = B_cache
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

    B_cache_0 <- if (logdet_method == "woodbury") {
      Q_mat <- if (!is.null(prior_0$Q_matrix)) prior_0$Q_matrix else
        if (is_matrix(prior_0$Q)) prior_0$Q else NULL
      CQ0 <- Matrix::Cholesky(Q_mat, LDL = FALSE, perm = TRUE)
      .make_B(CQ0)
    } else NULL

    phi_opt <- stats::optimize(
      f = function(lp) {
        res <- .eval_reml_ll(
          exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q0, prior_0,
          p, n, n_eff, logdet_Q0, probes, n_lanczos_steps,
          logdet_method, pcg_tol, pcg_maxit, x_warm,
          eval_constraint   = .eval_constraint,
          logdet_correction = .logdet_correction,
          solver            = solver,
          B_cache           = B_cache_0
        )
        if (!is.finite(res$ll)) Inf else -res$ll
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

  B_cache_opt <- if (logdet_method == "woodbury") {
    Q_mat <- if (!is.null(prior_opt$Q_matrix)) prior_opt$Q_matrix else
      if (is_matrix(prior_opt$Q)) prior_opt$Q else NULL
    CQ_opt <- Matrix::Cholesky(Q_mat, LDL = FALSE, perm = TRUE)
    .make_B(CQ_opt)
  } else NULL

  phi_final <- .profile_phi(
    AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, n_eff, logdet_Q_opt, probes, n_lanczos_steps,
    log_phi_lower, log_phi_upper,
    logdet_method, pcg_tol, pcg_maxit, x_warm,
    eval_constraint   = .eval_constraint,
    logdet_correction = .logdet_correction,
    solver            = solver,
    precond           = precond_apply,
    B_cache           = B_cache_opt
  )

  ll_final <- .eval_reml_ll(
    phi_final, AtRinvy, yRinvy, apply_AtRinvA, apply_Q_opt, prior_opt,
    p, n, n_eff, logdet_Q_opt, probes, n_lanczos_steps,
    logdet_method, pcg_tol_final, pcg_maxit, x_warm,
    eval_constraint   = .eval_constraint,
    logdet_correction = .logdet_correction,
    solver            = solver,
    precond           = precond_apply,
    B_cache           = B_cache_opt
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
                         eval_constraint   = NULL,
                         logdet_correction = NULL,
                         solver  = "pcg",
                         precond = NULL,
                         B_cache = NULL) {
  opt <- stats::optimize(
    f = function(lp) {
      res <- .eval_reml_ll(
        exp(lp), AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
        p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
        logdet_method, pcg_tol, pcg_maxit, x_warm,
        eval_constraint   = eval_constraint,
        logdet_correction = logdet_correction,
        solver            = solver,
        precond           = precond,
        B_cache           = B_cache
      )
      if (!is.finite(res$ll)) Inf else -res$ll
    },
    interval = c(log_phi_lower, log_phi_upper),
    tol      = 1e-4
  )
  exp(opt$minimum)
}

# -----------------------------------------------------------------------
# Internal: evaluate REML log likelihood at fixed (phi, theta)
# -----------------------------------------------------------------------
.eval_reml_ll <- function(phi, AtRinvy, yRinvy, apply_AtRinvA, apply_Q, prior,
                          p, n, n_eff, logdet_Q, probes, n_lanczos_steps,
                          logdet_method, pcg_tol, pcg_maxit, x_warm,
                          eval_constraint   = NULL,
                          logdet_correction = NULL,
                          solver  = "pcg",
                          precond = NULL,
                          B_cache = NULL) {

  apply_K <- function(v) apply_AtRinvA(v) + (1/phi) * apply_Q(v)

  # --- linear solve ---------------------------------------------------------
  if (solver == "cholesky" && logdet_method != "woodbury") {
    Q_mat <- if (!is.null(prior$Q_matrix)) prior$Q_matrix else
      if (is_matrix(prior$Q)) prior$Q else NULL
    if (!is.null(Q_mat) && !is.null(prior$AtRinvA_matrix)) {
      K_mat  <- Matrix::forceSymmetric(prior$AtRinvA_matrix + (1/phi) * Q_mat)
      chol_K <- tryCatch(
        Matrix::Cholesky(K_mat, LDL = FALSE, perm = TRUE),
        error = function(e) NULL
      )
      if (!is.null(chol_K)) {
        x <- as.numeric(Matrix::solve(chol_K, AtRinvy))
        return(.eval_reml_ll_from_x(
          x, phi, AtRinvy, yRinvy, apply_K, prior,
          p, n, n_eff, logdet_Q, probes, n_lanczos_steps, logdet_method,
          eval_constraint, logdet_correction,
          chol_K = chol_K, B_cache = B_cache
        ))
      }
    }
  }

  # PCG solve (projected operator when X_fixed present)
  pcg_res <- tryCatch(
    pcg(apply_K, AtRinvy, x0 = x_warm, precond = precond,
        tol = pcg_tol, maxit = pcg_maxit),
    error = function(e) list(converged = FALSE)
  )
  if (!pcg_res$converged) return(list(ll = -Inf, sigma2e = NA, x = x_warm))

  .eval_reml_ll_from_x(
    pcg_res$x, phi, AtRinvy, yRinvy, apply_K, prior,
    p, n, n_eff, logdet_Q, probes, n_lanczos_steps, logdet_method,
    eval_constraint, logdet_correction,
    chol_K = NULL, B_cache = B_cache
  )
}

# -----------------------------------------------------------------------
# Internal: complete ll given solution x
# -----------------------------------------------------------------------
.eval_reml_ll_from_x <- function(x, phi, AtRinvy, yRinvy, apply_K, prior,
                                 p, n, n_eff, logdet_Q, probes,
                                 n_lanczos_steps, logdet_method,
                                 eval_constraint,
                                 logdet_correction,
                                 chol_K  = NULL,
                                 B_cache = NULL) {

  # logdet_K: woodbury exact > cholesky exact > lanczos approximate
  logdet_K <- if (logdet_method == "woodbury" && !is.null(B_cache)) {
    .logdet_K_woodbury_inner(phi, B_cache, logdet_Q, p)
  } else if (!is.null(chol_K)) {
    as.numeric(Matrix::determinant(chol_K, logarithm = TRUE)$modulus) * 2
  } else {
    .eval_logdet_K(apply_K, prior, phi, p, probes, n_lanczos_steps,
                   logdet_method)
  }
  if (!is.finite(logdet_K)) return(list(ll = -Inf, sigma2e = NA, x = x))

  if (is.null(eval_constraint)) {
    # Quadratic form: y_proj'y_proj - (A'R^{-1}y_proj)' x
    # Already beta-profiled because AtRinvy uses projected y
    yHinvy <- yRinvy - as.numeric(crossprod(AtRinvy, x))
    if (yHinvy <= 0) return(list(ll = -Inf, sigma2e = NA, x = x))
    sigma2e <- yHinvy / n_eff

    ll <- -n_eff/2 * log(sigma2e) -
      1/2      * logdet_K     -
      p/2      * log(phi)     +
      1/2      * logdet_Q

    # Woodbury logdet already includes the full REML correction (exact).
    # Non-woodbury paths need the -1/2 log|X'H^{-1}X| correction separately.
    if (!is.null(logdet_correction)) {
      ldc <- logdet_correction(apply_K, p)
      if (!is.finite(ldc)) return(list(ll = -Inf, sigma2e = NA, x = x))
      ll <- ll + ldc
    }

  } else {
    constr <- eval_constraint(phi, environment(apply_K)$apply_Q, prior, n_eff)
    if (is.na(constr$sigma2e) || !is.finite(constr$logdet_M))
      return(list(ll = -Inf, sigma2e = NA, x = x))
    sigma2e <- constr$sigma2e
    ll <- -n_eff/2 * log(sigma2e) -
      1/2      * logdet_K     -
      p/2      * log(phi)     +
      1/2      * logdet_Q     -
      1/2      * constr$logdet_M
  }

  list(ll = ll, sigma2e = sigma2e, x = x)
}

# -----------------------------------------------------------------------
# Internal: woodbury logdet -- called per phi with cached B
# logdet(K_proj) = -p*log(phi) + logdet(Q) + logdet(I + phi*B)
# Note: logdet_Q is passed in separately; this returns logdet_K directly.
# -----------------------------------------------------------------------
.logdet_K_woodbury_inner <- function(phi, B, logdet_Q, p) {
  n   <- nrow(B)
  CM  <- tryCatch(chol(diag(n) + phi * B), error = function(e) NULL)
  if (is.null(CM)) return(-Inf)
  -p * log(phi) + logdet_Q + 2 * sum(log(diag(CM)))
}

# -----------------------------------------------------------------------
# Internal: log|K| via Lanczos or sparse Cholesky
# -----------------------------------------------------------------------
.eval_logdet_K <- function(apply_K, prior, phi, p,
                           probes, n_lanczos_steps, logdet_method) {
  if (logdet_method == "cholesky") {
    Q_mat <- if (!is.null(prior$Q_matrix)) prior$Q_matrix else
      if (is_matrix(prior$Q)) prior$Q else NULL
    if (!is.null(Q_mat) && !is.null(prior$AtRinvA_matrix)) {
      K_mat <- Matrix::forceSymmetric(prior$AtRinvA_matrix + (1/phi) * Q_mat)
      return(tryCatch(
        as.numeric(Matrix::determinant(K_mat, logarithm = TRUE)$modulus),
        error = function(e) lanczos_logdet(apply_K, probes, n_lanczos_steps)
      ))
    }
  }
  lanczos_logdet(apply_K, probes, n_lanczos_steps)
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
