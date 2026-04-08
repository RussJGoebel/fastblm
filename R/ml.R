#' Tune hyperparameters via marginal ML
#'
#' Maximizes the marginal ML likelihood over theta (shape parameters) and
#' profiles out phi (signal-to-noise) and sigma2e analytically at each theta.
#'
#' Unlike \code{\link{tune_reml}}, fixed effects are handled by including them
#' directly in the design matrix as an augmented system \eqn{[A | X_{\rm fixed}]}
#' with a flat (zero-precision) prior on the fixed-effect coefficients. This
#' avoids fixed-effect projection and keeps the system sparse, enabling exact
#' sparse Cholesky for both the linear solve and the log-determinant.
#'
#' This is ML not REML: the fixed-effect degrees of freedom are not subtracted
#' from n. With large n and small q this bias is negligible. See the package
#' vignette for a discussion of when REML vs ML matters.
#'
#' The log-likelihood at fixed (theta, phi) is:
#' \deqn{
#'   \ell = -\frac{n}{2}\log\sigma^2_e
#'          - \frac{1}{2}\log|K|
#'          - \frac{p}{2}\log\phi
#'          + \frac{1}{2}\log|Q|
#' }
#' where \eqn{K = A^\top R^{-1} A + \frac{1}{\phi} Q} (using the original
#' spatial A and Q, not the augmented versions) and
#' \eqn{\sigma^2_e = (y^\top R^{-1} y - \mu^\top A^\top R^{-1} y) / n}
#' with \eqn{\mu} the posterior mean of the augmented system.
#'
#' @param y numeric response vector length n
#' @param A n x p spatial design matrix (sparse dgCMatrix)
#' @param Q_fun function(theta) returning a list with elements:
#'   \code{Q} (p x p spatial prior precision, sparse matrix),
#'   \code{log_det_Q} (optional precomputed scalar logdet of Q).
#' @param X_fixed n x q matrix of fixed-effect covariates (optional).
#'   If supplied, columns are appended to A to form the augmented system.
#'   NULL means no fixed effects -- equivalent to \code{tune_reml} with
#'   Cholesky solver and logdet.
#' @param R_inv n x n inverse noise covariance. NULL = identity.
#' @param theta_init named numeric vector of initial theta values
#' @param lower lower bounds on theta
#' @param upper upper bounds on theta
#' @param log_phi_lower lower bound for phi search (log scale)
#' @param log_phi_upper upper bound for phi search (log scale)
#' @param verbose logical
#'
#' @return object of class \code{fastblm_tuned}, compatible with
#'   \code{\link{fit_fastblm}}, \code{\link{posterior_se}}, etc.
#' @export
tune_ml <- function(y, A, Q_fun,
                    X_fixed       = NULL,
                    R_inv         = NULL,
                    theta_init    = numeric(0),
                    lower         = rep(-Inf, length(theta_init)),
                    upper         = rep( Inf, length(theta_init)),
                    log_phi_lower = log(0.01),
                    log_phi_upper = log(1000),
                    verbose       = TRUE) {

  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(A)

  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")

  # -------------------------------------------------------------------------
  # Build augmented system if X_fixed supplied.
  # A_aug = [A | X_fixed],  Q_aug = bdiag(Q, 0_{q x q})
  # K_aug = A_aug' R^{-1} A_aug + (1/phi) Q_aug  -- sparse, Cholesky-able.
  #
  # The log-likelihood uses:
  #   logdet(K_aug) from the Cholesky factor (exact)
  #   logdet(Q_aug) = logdet(Q)  (zero block is improper, dropped by convention)
  #   -p_aug/2 * log(phi) has an extra -q/2*log(phi) vs the spatial-only term,
  #   but with q=2, p=20000 this is a shift of log(phi) << p/2*log(phi).
  #   We correct by using p not p_aug in the likelihood to match the spatial
  #   prior dimension.
  # -------------------------------------------------------------------------

  if (!is.null(X_fixed)) {
    X_fixed <- as.matrix(X_fixed)
    q       <- ncol(X_fixed)
    if (verbose)
      message(sprintf(
        "ML: augmenting A with %d fixed effect(s). n=%d, p=%d, p_aug=%d.",
        q, n, p, p + q))
    A_aug <- as(cbind(A, X_fixed), "dgCMatrix")
  } else {
    q     <- 0L
    A_aug <- A
  }

  p_aug <- ncol(A_aug)

  # --- precompute fixed quantities ------------------------------------------
  Rinv    <- resolve_Rinv(R_inv, n)
  Rinvy   <- as.numeric(Rinv %*% y)
  AtRinvy <- as.numeric(Matrix::crossprod(A_aug, Rinvy))   # p_aug vector
  yRinvy  <- as.numeric(crossprod(y, Rinvy))

  # Precompute A_aug' R^{-1} A_aug once -- reused across all theta/phi evals
  if (verbose) message("Forming A'R^{-1}A (once)...")
  AtRinvA_aug <- Matrix::crossprod(A_aug, Rinv %*% A_aug)   # p_aug x p_aug, sparse
  if (verbose)
    message(sprintf("  nnz=%d  density=%.3f%%",
                    nnzero(AtRinvA_aug),
                    nnzero(AtRinvA_aug) / p_aug^2 * 100))

  iter_count <- 0L
  history    <- list()

  # -------------------------------------------------------------------------
  # Core ll evaluation at fixed (theta, phi)
  # Returns list(ll, sigma2e, mu) where mu is the p_aug posterior mean.
  # -------------------------------------------------------------------------
  .eval_ml_ll <- function(phi, Q_aug, logdet_Q, chol_K = NULL) {
    # Form and factor K_aug if not pre-supplied
    if (is.null(chol_K)) {
      K_aug  <- Matrix::forceSymmetric(AtRinvA_aug + (1/phi) * Q_aug)
      chol_K <- tryCatch(
        Matrix::Cholesky(K_aug, LDL = FALSE, perm = TRUE),
        error = function(e) NULL
      )
      if (is.null(chol_K)) return(list(ll = -Inf, sigma2e = NA, mu = NULL))
    }

    mu      <- as.numeric(Matrix::solve(chol_K, AtRinvy))
    yHinvy  <- yRinvy - as.numeric(crossprod(AtRinvy, mu))
    if (yHinvy <= 0) return(list(ll = -Inf, sigma2e = NA, mu = mu))
    sigma2e <- yHinvy / n

    # logdet(K_aug) exact from Cholesky
    logdet_K <- as.numeric(
      Matrix::determinant(chol_K, logarithm = TRUE, sqrt = TRUE)$modulus) * 2

    # Use p (spatial dimension) not p_aug for the prior normalisation terms
    # -- the q extra beta columns have improper (flat) prior
    ll <- -n/2   * log(sigma2e) -
      1/2    * logdet_K     -
      p/2    * log(phi)     +
      1/2    * logdet_Q

    list(ll = ll, sigma2e = sigma2e, mu = mu, chol_K = chol_K)
  }

  # -------------------------------------------------------------------------
  # Profile phi at fixed theta via golden section
  # -------------------------------------------------------------------------
  .profile_phi_ml <- function(Q_aug, logdet_Q, x_warm) {
    opt <- stats::optimize(
      f = function(lp) {
        res <- .eval_ml_ll(exp(lp), Q_aug, logdet_Q)
        if (!is.finite(res$ll)) Inf else -res$ll
      },
      interval = c(log_phi_lower, log_phi_upper),
      tol      = 1e-4
    )
    exp(opt$minimum)
  }

  # -------------------------------------------------------------------------
  # Build Q_aug from Q_fun output
  # -------------------------------------------------------------------------
  .make_Q_aug <- function(prior) {
    Q_sp <- prior$Q
    if (!is_matrix(Q_sp))
      stop("tune_ml requires Q_fun to return an explicit sparse matrix for Q.")

    if (q > 0L) {
      Q_aug <- Matrix::forceSymmetric(Matrix::bdiag(
        Q_sp,
        Matrix::Matrix(0, nrow = q, ncol = q)
      ))
    } else {
      Q_aug <- Q_sp
    }

    # logdet of Q: spatial block only (zero block is improper, convention = drop)
    logdet_Q <- if (!is.null(prior$log_det_Q)) {
      prior$log_det_Q
    } else {
      CQ <- tryCatch(
        Matrix::Cholesky(Q_sp, LDL = FALSE, perm = TRUE),
        error = function(e) NULL
      )
      if (is.null(CQ)) return(NULL)
      as.numeric(
        Matrix::determinant(CQ, logarithm = TRUE, sqrt = TRUE)$modulus) * 2
    }

    list(Q_aug = Q_aug, logdet_Q = logdet_Q)
  }

  # -------------------------------------------------------------------------
  # Outer objective over theta
  # -------------------------------------------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)

    qa <- .make_Q_aug(prior)
    if (is.null(qa)) return(.Machine$double.xmax)

    phi_hat <- .profile_phi_ml(qa$Q_aug, qa$logdet_Q, x_warm = NULL)
    ll_res  <- .eval_ml_ll(phi_hat, qa$Q_aug, qa$logdet_Q)
    if (!is.finite(ll_res$ll)) return(.Machine$double.xmax)

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

  # -------------------------------------------------------------------------
  # Dispatch optimizer
  # -------------------------------------------------------------------------
  if (verbose) message("Starting ML optimisation...")

  if (length(theta_init) == 0L) {
    prior_0 <- Q_fun(numeric(0))
    qa_0    <- .make_Q_aug(prior_0)
    phi_opt <- stats::optimize(
      f        = function(lp) {
        res <- .eval_ml_ll(exp(lp), qa_0$Q_aug, qa_0$logdet_Q)
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

  # -------------------------------------------------------------------------
  # Recover at optimum
  # -------------------------------------------------------------------------
  prior_opt <- Q_fun(theta_opt)
  qa_opt    <- .make_Q_aug(prior_opt)

  phi_final <- .profile_phi_ml(qa_opt$Q_aug, qa_opt$logdet_Q, x_warm = NULL)
  ll_final  <- .eval_ml_ll(phi_final, qa_opt$Q_aug, qa_opt$logdet_Q)

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
      Q       = prior_opt$Q,           # spatial Q only (not augmented)
      value   = ll_final$ll,
      method  = "ml",
      optim   = optim_res,
      history = hist_df
    ),
    class = "fastblm_tuned"
  )
}
