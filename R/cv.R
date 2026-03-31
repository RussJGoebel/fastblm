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
#' @param pcg_tol PCG tolerance
#' @param pcg_maxit max PCG iterations
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
                    pcg_tol       = 1e-6,
                    pcg_maxit     = NULL,
                    verbose       = TRUE) {

  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(A)
  if (is.null(pcg_maxit)) pcg_maxit <- 4L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")

  score_fn <- .make_score_fn(score)

  # --- make folds -------------------------------------------------------
  folds <- .make_folds(n, k)

  # --- iteration history ------------------------------------------------
  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta ---------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- tryCatch(Q_fun(theta), error = function(e) NULL)
    if (is.null(prior) || is.null(prior$Q)) return(.Machine$double.xmax)

    # profile phi via CV
    phi_hat <- .profile_phi_cv(
      y, A, prior, folds, score_fn,
      log_phi_lower, log_phi_upper,
      pcg_tol, pcg_maxit
    )

    cv_score <- .eval_cv(
      y, A, prior, phi_hat, folds, score_fn,
      pcg_tol, pcg_maxit
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

  # --- dispatch optimizer -----------------------------------------------
  if (verbose) message("Starting CV optimisation...")

  if (length(theta_init) == 0L) {
    prior_0  <- Q_fun(numeric(0))
    phi_opt  <- exp(stats::optimize(
      f        = function(lp) .eval_cv(y, A, prior_0, exp(lp), folds,
                                       score_fn, pcg_tol, pcg_maxit),
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

  # --- recover at optimum -----------------------------------------------
  prior_opt <- Q_fun(theta_opt)
  phi_final <- .profile_phi_cv(
    y, A, prior_opt, folds, score_fn,
    log_phi_lower, log_phi_upper,
    pcg_tol, pcg_maxit
  )

  # fit on full data at optimum to get sigma2e
  fit_full <- fit_fastblm(y, A, prior_opt$Q, phi = phi_final,
                          solver = "pcg",
                          pcg_tol = pcg_tol, pcg_maxit = pcg_maxit)

  cv_final <- .eval_cv(y, A, prior_opt, phi_final, folds, score_fn,
                       pcg_tol, pcg_maxit)

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
                            pcg_tol, pcg_maxit) {
  opt <- stats::optimize(
    f        = function(lp) .eval_cv(y, A, prior, exp(lp), folds,
                                     score_fn, pcg_tol, pcg_maxit),
    interval = c(log_phi_lower, log_phi_upper),
    tol      = 1e-4
  )
  exp(opt$minimum)
}

# -----------------------------------------------------------------------
# Internal: evaluate k-fold CV score at fixed (theta, phi)
# -----------------------------------------------------------------------
.eval_cv <- function(y, A, prior, phi, folds, score_fn, pcg_tol, pcg_maxit) {
  k      <- max(folds)
  scores <- numeric(k)

  for (fold in seq_len(k)) {
    test_idx  <- which(folds == fold)
    train_idx <- which(folds != fold)

    y_train <- y[train_idx]
    A_train <- A[train_idx, , drop = FALSE]
    y_test  <- y[test_idx]
    A_test  <- A[test_idx, , drop = FALSE]

    fit_fold <- tryCatch(
      fit_fastblm(y_train, A_train, prior$Q, phi = phi,
                  solver = "pcg", pcg_tol = pcg_tol, pcg_maxit = pcg_maxit),
      error = function(e) NULL
    )
    if (is.null(fit_fold)) return(Inf)

    y_pred    <- as.numeric(A_test %*% fit_fold$posterior_mean)
    scores[fold] <- score_fn(y_test, y_pred)
  }

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
