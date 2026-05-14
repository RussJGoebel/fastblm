#' Tune hyperparameters via k-fold cross validation
#'
#' Selects theta and phi by minimizing k-fold CV prediction error.
#' At each candidate (theta, phi), fits on training folds and evaluates
#' prediction error on held-out fold.
#'
#' @param y numeric response vector length n
#' @param A n x p design matrix
#' @param Q_fun function returning prior. Two calling conventions are supported:
#'   \itemize{
#'     \item \code{function(theta)}: Q does not depend on the fold. Called once
#'       per theta evaluation and reused across all folds and phi evaluations.
#'     \item \code{function(theta, A_train)}: Q depends on the training data
#'       (e.g. soft RSR where \eqn{C = X_{\rm obs,train}^\top A_{\rm train}}
#'       is fold-specific). Called once per fold per phi evaluation.
#'   }
#'   In both cases the return value must be a list with element \code{Q}
#'   (required) and optionally \code{log_det_Q}.
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
#' @param constraint optional constraint specification. One of:
#'   \itemize{
#'     \item \code{NULL} (default): no constraint applied.
#'     \item a \code{q x p} matrix \code{C}: the fixed constraint \eqn{C\gamma = 0}
#'       applied identically on every fold.
#'     \item a function \code{function(train_idx, A_train) -> C}: called once per
#'       fold to produce a fold-specific constraint matrix.
#'   }
#' @param precond_fun optional function \code{function(phi, prior, A_train, y_train)}
#'   returning a preconditioner function for PCG. NULL = identity.
#' @param optim_method optimizer to use when \code{length(theta_init) > 1}.
#'   One of:
#'   \itemize{
#'     \item \code{"lbfgsb"} (default): L-BFGS-B via \code{stats::optim}.
#'     \item \code{"coordinate"}: Coordinate-wise golden section search.
#'   }
#' @param coord_tol convergence tolerance for coordinate descent. Default \code{1e-3}.
#' @param coord_maxit maximum number of coordinate descent cycles. Default 20.
#' @param lbfgsb_control named list of control parameters passed to
#'   \code{stats::optim} when \code{optim_method = "lbfgsb"}.
#' @param weighted_folds logical. If TRUE, CV scores are weighted by fold size
#'   (number of observations per fold) rather than averaging folds equally.
#'   Useful when folds are unequal in size, e.g. with spatially blocked CV.
#'   Default FALSE preserves existing behaviour.
#' @param seed optional integer seed for fold assignment reproducibility.
#' @param folds optional integer vector of length n giving pre-specified fold
#'   assignments. If supplied, \code{k} and \code{seed} are ignored.
#' @param parallel logical. If TRUE, evaluates CV folds in parallel. Default FALSE.
#' @param verbose logical
#'
#' @return object of class fastblm_tuned
#' @export
tune_cv <- function(y, A, Q_fun,
                    R_inv          = NULL,
                    theta_init     = numeric(0),
                    lower          = rep(-Inf, length(theta_init)),
                    upper          = rep( Inf, length(theta_init)),
                    log_phi_lower  = log(0.01),
                    log_phi_upper  = log(1000),
                    k              = 5L,
                    score          = "mse",
                    solver         = c("cholesky", "pcg", "woodbury"),
                    pcg_tol        = 1e-6,
                    pcg_maxit      = NULL,
                    Q_inv          = NULL,
                    constraint     = NULL,
                    precond_fun    = NULL,
                    optim_method   = c("lbfgsb", "coordinate"),
                    coord_tol      = 1e-3,
                    coord_maxit    = 20L,
                    lbfgsb_control = list(),
                    weighted_folds = FALSE,
                    seed           = NULL,
                    folds          = NULL,
                    parallel       = FALSE,
                    verbose        = TRUE) {

  y            <- as.numeric(y)
  n            <- length(y)
  p            <- ncol(A)
  solver       <- match.arg(solver)
  optim_method <- match.arg(optim_method)

  if (parallel && !requireNamespace("future.apply", quietly = TRUE))
    stop("parallel = TRUE requires the 'future.apply' package.")

  if (is.null(pcg_maxit)) pcg_maxit <- 2L * p
  if (length(theta_init) > 0L && is.null(names(theta_init)))
    stop("`theta_init` must be a named vector.")
  if (solver == "woodbury" && is.null(Q_inv))
    stop("solver = 'woodbury' requires Q_inv.")

  score_fn <- .make_score_fn(score)

  constraint_fn <- .make_constraint_fn(constraint, p)

  # Detect whether Q_fun is fold-aware (2 args) or global (1 arg).
  q_fun_fold_aware <- (length(formals(Q_fun)) >= 2L)

  # --- make or validate folds ----------------------------------------------
  if (!is.null(folds)) {
    if (length(folds) != n) stop("`folds` must have length n.")
    if (anyNA(folds))       stop("`folds` contains NA values.")
    folds <- as.integer(folds)
    k     <- max(folds)
  } else {
    if (!is.null(seed)) set.seed(seed)
    folds <- .make_folds(n, k)
  }

  # --- fold weights --------------------------------------------------------
  # Proportional to fold size, normalised so weights sum to k.
  # When all folds are equal size, weights are all 1 and weighted.mean = mean.
  fold_weights <- if (weighted_folds) {
    w <- tabulate(folds)[seq_len(k)]
    w * k / sum(w)
  } else {
    NULL
  }

  if (weighted_folds && verbose) {
    message(sprintf("weighted_folds = TRUE: fold sizes = %s",
                    paste(tabulate(folds)[seq_len(k)], collapse = " ")))
  }

  # --- iteration history ---------------------------------------------------
  iter_count <- 0L
  history    <- list()

  # --- outer objective over theta ------------------------------------------
  .objective <- function(theta) {
    if (length(theta) > 0L) names(theta) <- names(theta_init)
    iter_count <<- iter_count + 1L

    prior <- if (!q_fun_fold_aware) {
      tryCatch(Q_fun(theta), error = function(e) NULL)
    } else NULL

    if (!q_fun_fold_aware && (is.null(prior) || is.null(prior$Q)))
      return(.Machine$double.xmax)

    fold_C_list <- .precompute_fold_constraints(constraint_fn, A, folds)

    phi_hat <- .profile_phi_cv(
      y, A, Q_fun, theta, prior, folds, score_fn,
      log_phi_lower, log_phi_upper,
      solver, pcg_tol, pcg_maxit, Q_inv,
      R_inv        = R_inv,
      fold_C_list  = fold_C_list,
      precond_fun  = precond_fun,
      fold_weights = fold_weights,
      parallel     = parallel
    )

    cv_score <- .eval_cv(
      y, A, Q_fun, theta, prior, phi_hat, folds, score_fn,
      solver, pcg_tol, pcg_maxit, Q_inv,
      R_inv        = R_inv,
      fold_C_list  = fold_C_list,
      precond_fun  = precond_fun,
      fold_weights = fold_weights,
      parallel     = parallel
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
    prior_0     <- if (!q_fun_fold_aware) Q_fun(numeric(0)) else NULL
    fold_C_list <- .precompute_fold_constraints(constraint_fn, A, folds)
    phi_opt <- exp(stats::optimize(
      f        = function(lp) .eval_cv(
        y, A, Q_fun, numeric(0), prior_0, exp(lp), folds,
        score_fn, solver, pcg_tol, pcg_maxit, Q_inv,
        R_inv        = R_inv,
        fold_C_list  = fold_C_list,
        precond_fun  = precond_fun,
        fold_weights = fold_weights,
        parallel     = parallel),
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

  } else if (optim_method == "coordinate") {
    theta_cur <- unname(theta_init)
    names(theta_cur) <- names(theta_init)

    optim_res <- list(convergence = 1L)
    for (cycle in seq_len(coord_maxit)) {
      theta_prev <- theta_cur
      for (j in seq_along(theta_cur)) {
        opt_j <- stats::optimize(
          f        = function(val) {
            theta_try <- theta_cur
            theta_try[j] <- val
            names(theta_try) <- names(theta_init)
            .objective(theta_try)
          },
          interval = c(unname(lower)[j], unname(upper)[j]),
          tol      = 1e-4
        )
        theta_cur[j] <- opt_j$minimum
      }
      max_change <- max(abs(theta_cur - theta_prev))
      if (verbose)
        message(sprintf("  [coordinate cycle %d]  max_change=%.4g", cycle, max_change))
      if (max_change < coord_tol) {
        optim_res <- list(convergence = 0L, cycles = cycle)
        break
      }
    }
    theta_opt <- theta_cur
    names(theta_opt) <- names(theta_init)

  } else {
    default_control <- list(maxit = 300L, trace = if (verbose) 1L else 0L)
    ctl <- modifyList(default_control, lbfgsb_control)
    if (!is.null(ctl$ndeps))
      ctl$ndeps <- rep_len(ctl$ndeps, length(theta_init))
    opt <- stats::optim(
      par     = unname(theta_init),
      fn      = function(par) {
        theta <- par; names(theta) <- names(theta_init)
        .objective(theta)
      },
      method  = "L-BFGS-B",
      lower   = unname(lower),
      upper   = unname(upper),
      control = ctl
    )
    theta_opt        <- opt$par
    names(theta_opt) <- names(theta_init)
    optim_res        <- opt
  }

  # --- recover at optimum --------------------------------------------------
  if (length(theta_opt) > 0L) {
    prior_opt   <- if (!q_fun_fold_aware) Q_fun(theta_opt) else NULL
    fold_C_list <- .precompute_fold_constraints(constraint_fn, A, folds)
    phi_final   <- .profile_phi_cv(
      y, A, Q_fun, theta_opt, prior_opt, folds, score_fn,
      log_phi_lower, log_phi_upper,
      solver, pcg_tol, pcg_maxit, Q_inv,
      R_inv        = R_inv,
      fold_C_list  = fold_C_list,
      precond_fun  = precond_fun,
      fold_weights = fold_weights,
      parallel     = parallel
    )
  } else {
    prior_opt <- if (!q_fun_fold_aware) Q_fun(numeric(0)) else NULL
    phi_final <- phi_opt
  }

  prior_full <- if (q_fun_fold_aware) Q_fun(theta_opt, A) else prior_opt

  fit_full <- fit_fastblm(y, A, prior_full$Q, phi = phi_final,
                          R_inv       = R_inv,
                          solver      = solver,
                          Q_inv       = Q_inv,
                          pcg_tol     = pcg_tol,
                          pcg_maxit   = pcg_maxit,
                          pcg_precond = prior_full$precond)

  cv_final <- .eval_cv(y, A, Q_fun, theta_opt, prior_opt, phi_final, folds,
                       score_fn, solver, pcg_tol, pcg_maxit, Q_inv,
                       R_inv        = R_inv,
                       fold_C_list  = fold_C_list,
                       precond_fun  = precond_fun,
                       fold_weights = fold_weights,
                       parallel     = parallel)

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
      Q       = prior_full$Q,
      value   = cv_final,
      method  = paste0("cv_", k, "fold_", score),
      optim   = optim_res,
      history = hist_df
    ),
    class = "fastblm_tuned"
  )
}

# -----------------------------------------------------------------------
# Internal: normalise constraint to a function or NULL
# -----------------------------------------------------------------------
.make_constraint_fn <- function(constraint, p) {
  if (is.null(constraint)) return(NULL)
  if (is.function(constraint)) {
    if (length(formals(constraint)) != 2L)
      stop("`constraint` function must have exactly two arguments: function(train_idx, A_train).")
    return(constraint)
  }
  if (is_matrix(constraint)) {
    if (ncol(constraint) != p)
      stop(sprintf("`constraint` matrix must have ncol = p = %d.", p))
    C <- constraint
    return(function(train_idx, A_train) C)
  }
  stop("`constraint` must be NULL, a matrix, or a function(train_idx, A_train).")
}

# -----------------------------------------------------------------------
# Internal: precompute per-fold constraint matrices
# -----------------------------------------------------------------------
.precompute_fold_constraints <- function(constraint_fn, A, folds) {
  if (is.null(constraint_fn)) return(NULL)
  k <- max(folds)
  lapply(seq_len(k), function(fold) {
    train_idx <- which(folds != fold)
    constraint_fn(train_idx, A[train_idx, , drop = FALSE])
  })
}

# -----------------------------------------------------------------------
# Internal: profile phi via CV at fixed theta
# -----------------------------------------------------------------------
.profile_phi_cv <- function(y, A, Q_fun, theta, prior, folds, score_fn,
                            log_phi_lower, log_phi_upper,
                            solver, pcg_tol, pcg_maxit, Q_inv,
                            R_inv        = NULL,
                            fold_C_list  = NULL,
                            precond_fun  = NULL,
                            fold_weights = NULL,
                            parallel     = FALSE) {
  opt <- stats::optimize(
    f        = function(lp) .eval_cv(
      y, A, Q_fun, theta, prior, exp(lp), folds,
      score_fn, solver, pcg_tol, pcg_maxit, Q_inv,
      R_inv        = R_inv,
      fold_C_list  = fold_C_list,
      precond_fun  = precond_fun,
      fold_weights = fold_weights,
      parallel     = parallel),
    interval = c(log_phi_lower, log_phi_upper),
    tol      = 1e-4
  )
  exp(opt$minimum)
}

# -----------------------------------------------------------------------
# Internal: evaluate k-fold CV score at fixed (theta, phi)
# -----------------------------------------------------------------------
.eval_cv <- function(y, A, Q_fun, theta, prior, phi, folds, score_fn,
                     solver, pcg_tol, pcg_maxit, Q_inv,
                     R_inv        = NULL,
                     fold_C_list  = NULL,
                     precond_fun  = NULL,
                     fold_weights = NULL,
                     parallel     = FALSE) {
  k                <- max(folds)
  q_fun_fold_aware <- (length(formals(Q_fun)) >= 2L)

  .fit_fold <- function(fold) {
    test_idx  <- which(folds == fold)
    train_idx <- which(folds != fold)

    y_train <- y[train_idx]
    A_train <- A[train_idx, , drop = FALSE]
    y_test  <- y[test_idx]
    A_test  <- A[test_idx, , drop = FALSE]

    R_inv_train <- if (!is.null(R_inv)) R_inv[train_idx, train_idx] else NULL

    fold_prior <- if (q_fun_fold_aware)
      tryCatch(Q_fun(theta, A_train), error = function(e) NULL)
    else
      prior

    if (is.null(fold_prior) || is.null(fold_prior$Q)) return(Inf)

    pcg_precond <- if (!is.null(precond_fun))
      precond_fun(phi, fold_prior, A_train, y_train)
    else
      fold_prior$precond

    fit_fold <- tryCatch(
      fit_fastblm(y_train, A_train, fold_prior$Q, phi = phi,
                  R_inv       = R_inv_train,
                  solver      = solver,
                  Q_inv       = Q_inv,
                  pcg_tol     = pcg_tol,
                  pcg_maxit   = pcg_maxit,
                  pcg_precond = pcg_precond),
      error = function(e) {
        warning(sprintf("fold %d: fit_fastblm failed at phi=%.4g: %s",
                        fold, phi, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit_fold)) return(Inf)

    if (!is.null(fold_C_list)) {
      fit_fold <- tryCatch(
        constrain(fit_fold, fold_C_list[[fold]]),
        error = function(e) {
          warning(sprintf("fold %d: constrain failed at phi=%.4g: %s",
                          fold, phi, conditionMessage(e)))
          NULL
        }
      )
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
        y                = y,
        A                = A,
        Q_fun            = Q_fun,
        theta            = theta,
        prior            = prior,
        folds            = folds,
        phi              = phi,
        R_inv            = R_inv,
        score_fn         = score_fn,
        solver           = solver,
        Q_inv            = Q_inv,
        pcg_tol          = pcg_tol,
        pcg_maxit        = pcg_maxit,
        fold_C_list      = fold_C_list,
        precond_fun      = precond_fun,
        q_fun_fold_aware = q_fun_fold_aware,
        fit_fastblm      = fastblm::fit_fastblm,
        constrain        = fastblm::constrain
      ),
      future.packages = c("Matrix", "fastblm")
    ))
  } else {
    vapply(seq_len(k), .fit_fold, numeric(1L))
  }

  if (any(!is.finite(scores))) return(Inf)

  if (!is.null(fold_weights)) {
    weighted.mean(scores, fold_weights)
  } else {
    mean(scores)
  }
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
