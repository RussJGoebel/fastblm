#' Tune hyperparameters via k-fold cross validation
#'
#' Selects theta and phi by minimizing k-fold CV prediction error.
#' At each candidate (theta, phi), fits on training folds and evaluates
#' prediction error on held-out fold.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param Q_fun function(theta) -> list with Q (required), log_det_Q (optional)
#' @param R_inv n x n inverse noise covariance. NULL = identity.
#' @param theta_init named numeric vector of initial theta values
#' @param lower lower bounds on theta
#' @param upper upper bounds on theta
#' @param log_phi_lower lower bound for phi search (log scale)
#' @param log_phi_upper upper bound for phi search (log scale)
#' @param k number of folds. Default 5.
#' @param score scoring function. Default "mse". Options: "mse", "mae"
#' @param solver solver to use for each fold fit. One of \code{"cholesky"},
#'   \code{"pcg"}, or \code{"woodbury"}. Default \code{"cholesky"}.
#' @param pcg_tol PCG tolerance (only used when solver = "pcg")
#' @param pcg_maxit max PCG iterations (only used when solver = "pcg")
#' @param Q_inv optional p x p prior covariance or function v -> Q^{-1}v.
#'   Required when solver = "woodbury".
#' @param constraint_matrix optional constraint matrix C such that C \%*\% gamma = 0
#'   is enforced on the posterior. Supply either this or \code{X_cov}, not both.
#' @param X_cov optional n x q matrix of covariates used to compute the constraint
#'   matrix per fold as \code{t(X_cov[train,]) \%*\% A_train}. Supply either this
#'   or \code{constraint_matrix}, not both.
#' @param seed optional integer seed for fold assignment reproducibility. Ignored
#'   if \code{folds} is supplied.
#' @param folds optional integer vector of length n giving pre-specified fold
#'   assignments. If supplied, \code{k} and \code{seed} are ignored.
#' @param parallel logical. If TRUE, evaluates CV folds in parallel using
#'   \code{future.apply::future_lapply}. The parallelization backend is
#'   controlled by the user via \code{future::plan()} before calling this
#'   function. Default FALSE.
#' @param verbose logical
#'
#' @return object of class fastblm_tuned
#' @export
tune_cv <- function(y, A, Q_fun,
                    R_inv         = NULL,
                    theta_init    = numeric(0),
                    lower         = rep(-Inf, length(theta_init)),
                    upper         = rep( Inf, length(theta_init)),
                    log_phi_lower = log(0.01),
                    log_phi_upper = log(1000),
                    k             = 5L,
                    score         = "mse",
                    solver        = c("cholesky", "pcg", "woodbury"),
                    pcg_tol       = 1e-6,
                    pcg_maxit     = NULL,
                    Q_inv         = NULL,
                    constraint_matrix = NULL,
                    X_cov         = NULL,
                    seed          = NULL,
                    folds         = NULL,
                    parallel      = FALSE,
                    verbose       = TRUE) {

  y      <- as.numeric(y)
  n      <- length(y)
  p      <- ncol(A)
  solver <- match.arg(solver)

  if (parallel && !requireNamespace("future.apply", quietly = TRUE))
    stop("parallel = TRUE requires the 'future.apply' package. Install it with install.packages('future.apply').")

  if (is.null(pcg_maxit)) pcg_maxit <- 4L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")
  if (solver == "woodbury" && is.null(Q_inv))
    stop("solver = 'woodbury' requires Q_inv.")

  score_fn <- .make_score_fn(score)

  # --- validate constraint -------------------------------------------------
  use_constraint <- !is.null(constraint_matrix) || !is.null(X_cov)
  if (!is.null(constraint_matrix) && !is.null(X_cov))
    stop("Supply either `constraint_matrix` or `X_cov`, not both.")
  if (!is.null(constraint_matrix) && !is_matrix(constraint_matrix))
    stop("`constraint_matrix` must be a matrix.")
  if (!is.null(X_cov) && nrow(X_cov) != n)
    stop("`X_cov` must have nrow = n.")

  # --- make or validate folds ----------------------------------------------
  if (!is.null(folds)) {
    if (length(folds) != n)
      stop("`folds` must have length n.")
    if (anyNA(folds))
      stop("`folds` contains NA values.")
    folds <- as.integer(folds)
    k     <- max(folds)
  } else {
    if (!is.null(seed)) set.seed(seed)
    folds <- .make_folds(n, k)
  }

  # --- iteration history ---------------------------------------------------
  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta ------------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)

    phi_hat <- .profile_phi_cv(
      y, A, prior, folds, score_fn,
      log_phi_lower, log_phi_upper,
      solver, pcg_tol, pcg_maxit, Q_inv,
      constraint_matrix = constraint_matrix, X_cov = X_cov,
      parallel = parallel
    )

    cv_score <- .eval_cv(
      y, A, prior, phi_hat, folds, score_fn,
      solver, pcg_tol, pcg_maxit, Q_inv,
      constraint_matrix = constraint_matrix, X_cov = X_cov,
      parallel = parallel
    )

    if (!is.finite(cv_score)) return(.Machine$double.xmax)

    if (verbose) {
      theta_str <- if (length(theta) > 0L)
        paste(sprintf("%s=%.4g", names(theta), theta), collapse = "  ")
      else ""
      message(sprintf("  iter %d: %s  log_phi=%.3f  cv_%s=%.4g",
                      iter_count, theta_str, log(phi_hat), score, cv_score))
    }

    history[[iter_count]] <<- list(
      theta    = theta,
      phi      = phi_hat,
      cv_score = cv_score
    )

    cv_score
  }

  # --- dispatch optimizer --------------------------------------------------
  if (verbose) message("Starting CV optimisation...")

  if (length(theta_init) == 0L) {
    prior_0 <- Q_fun(numeric(0))
    phi_opt <- exp(stats::optimize(
      f        = function(lp) .eval_cv(y, A, prior_0, exp(lp), folds,
                                       score_fn, solver, pcg_tol, pcg_maxit, Q_inv,
                                       constraint_matrix = constraint_matrix,
                                       X_cov = X_cov,
                                       parallel = parallel),
      interval = c(log_phi_lower, log_phi_upper),
      tol      = 1e-4
    )$minimum)
    theta_opt <- numeric(0)
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

  # --- recover at optimum --------------------------------------------------
  prior_opt <- Q_fun(theta_opt)
  phi_final <- .profile_phi_cv(
    y, A, prior_opt, folds, score_fn,
    log_phi_lower, log_phi_upper,
    solver, pcg_tol, pcg_maxit, Q_inv,
    constraint_matrix = constraint_matrix, X_cov = X_cov,
    parallel = parallel
  )

  fit_full <- fit_fastblm(y, A, prior_opt$Q, phi = phi_final,
                          solver    = solver,
                          Q_inv     = Q_inv,
                          pcg_tol   = pcg_tol,
                          pcg_maxit = pcg_maxit)

  cv_final <- .eval_cv(y, A, prior_opt, phi_final, folds, score_fn,
                       solver, pcg_tol, pcg_maxit, Q_inv,
                       constraint_matrix = constraint_matrix, X_cov = X_cov,
                       parallel = parallel)

  if (verbose)
    message(sprintf("Optimum: phi=%.4g  sigma2e=%.4g  cv_%s=%.4g",
                    phi_final, fit_full$sigma2e, score, cv_final))

  hist_df <- if (length(history) > 0L)
    do.call(rbind, lapply(history, as.data.frame))
  else data.frame()

  structure(
    list(
      theta   = theta_opt,
      phi     = phi_final,
      sigma2e = fit_full$sigma2e,
      sigma2b = phi_final * fit_full$sigma2e,
      Q       = prior_opt$Q,
      value   = cv_final,
      method  = paste0("cv_", k, "fold_", score),
      optim   = optim_res,
      history = hist_df
    ),
    class = "fastblm_tuned"
  )
}

# -----------------------------------------------------------------------
# Internal: profile phi via CV at fixed theta
# -----------------------------------------------------------------------
.profile_phi_cv <- function(y, A, prior, folds, score_fn,
                            log_phi_lower, log_phi_upper,
                            solver, pcg_tol, pcg_maxit, Q_inv,
                            constraint_matrix = NULL, X_cov = NULL,
                            parallel = FALSE) {
  opt <- stats::optimize(
    f        = function(lp) .eval_cv(y, A, prior, exp(lp), folds,
                                     score_fn, solver, pcg_tol, pcg_maxit, Q_inv,
                                     constraint_matrix = constraint_matrix,
                                     X_cov = X_cov,
                                     parallel = parallel),
    interval = c(log_phi_lower, log_phi_upper),
    tol      = 1e-4
  )
  exp(opt$minimum)
}

# -----------------------------------------------------------------------
# Internal: evaluate k-fold CV score at fixed (theta, phi)
# -----------------------------------------------------------------------
.eval_cv <- function(y, A, prior, phi, folds, score_fn,
                     solver, pcg_tol, pcg_maxit, Q_inv,
                     constraint_matrix = NULL, X_cov = NULL,
                     parallel = FALSE) {
  k <- max(folds)

  .fit_fold <- function(fold) {
    test_idx  <- which(folds == fold)
    train_idx <- which(folds != fold)

    y_train <- y[train_idx]
    A_train <- A[train_idx, , drop = FALSE]
    y_test  <- y[test_idx]
    A_test  <- A[test_idx, , drop = FALSE]

    fit_fold <- tryCatch(
      fit_fastblm(y_train, A_train, prior$Q, phi = phi,
                  solver    = solver,
                  Q_inv     = Q_inv,
                  pcg_tol   = pcg_tol,
                  pcg_maxit = pcg_maxit),
      error = function(e) NULL
    )
    if (is.null(fit_fold)) return(Inf)

    if (!is.null(X_cov)) {
      C_fold   <- t(X_cov[train_idx, , drop = FALSE]) %*% A_train
      fit_fold <- tryCatch(constrain(fit_fold, C_fold), error = function(e) NULL)
      if (is.null(fit_fold)) return(Inf)
    } else if (!is.null(constraint_matrix)) {
      fit_fold <- tryCatch(constrain(fit_fold, constraint_matrix),
                           error = function(e) NULL)
      if (is.null(fit_fold)) return(Inf)
    }

    y_pred <- as.numeric(A_test %*% fit_fold$posterior_mean)
    score_fn(y_test, y_pred)
  }

  scores <- if (parallel) {
    unlist(future.apply::future_lapply(
      seq_len(k),
      .fit_fold,
      future.seed    = TRUE,
      future.globals = list(
        y                 = y,
        A                 = A,
        folds             = folds,
        prior             = prior,
        phi               = phi,
        score_fn          = score_fn,
        solver            = solver,
        Q_inv             = Q_inv,
        pcg_tol           = pcg_tol,
        pcg_maxit         = pcg_maxit,
        X_cov             = X_cov,
        constraint_matrix = constraint_matrix,
        fit_fastblm       = fastblm::fit_fastblm,
        constrain         = fastblm::constrain
      ),
      future.packages = c("Matrix", "fastblm")
    ))
  } else {
    vapply(seq_len(k), .fit_fold, numeric(1L))
  }

  if (any(!is.finite(scores))) return(Inf)
  mean(scores)
}

# -----------------------------------------------------------------------
# Internal: make fold assignments
# -----------------------------------------------------------------------
.make_folds <- function(n, k) {
  folds <- rep_len(seq_len(k), n)
  sample(folds)
}

# -----------------------------------------------------------------------
# Internal: score functions
# -----------------------------------------------------------------------
.make_score_fn <- function(score) {
  switch(score,
         mse = function(y, yhat) mean((y - yhat)^2),
         mae = function(y, yhat) mean(abs(y - yhat)),
         stop(sprintf("Unknown score '%s'. Options: 'mse', 'mae'.", score))
  )
}
