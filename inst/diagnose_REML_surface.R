# =============================================================================
# diagnose_reml_surface.R
#
# Plots the REML log-likelihood surface over a grid of (rho, log_phi) values
# for a single simulated dataset. Helps diagnose whether phi=999 blowup is
# a flat ridge, a code issue, or an optimizer problem.
#
# Run AFTER sourcing utils.R, solvers.R, fit.R, reml.R, and after
# test_rho_recovery.R has been sourced (reuses its helpers).
# =============================================================================

# ---- reproduce the A_direct rho=0.3 dataset ---------------------------------
set.seed(42)

m       <- 20
p       <- m * m
phi     <- 5.0
sigma2e <- 1.0
W       <- make_grid_W(m)

true_rho <- 0.3
x        <- simulate_sar(W, true_rho, phi, sigma2e)
y        <- x + rnorm(p, sd = sqrt(sigma2e))
A        <- as.matrix(Matrix::Diagonal(p))

# ---- precompute once --------------------------------------------------------
Q_fun   <- make_sar_Q_fun(W)
n       <- length(y)
Rinvy   <- y                                         # R = I
AtRinvy <- as.numeric(crossprod(A, y))
yRinvy  <- as.numeric(crossprod(y))
AtRinvA <- Matrix::crossprod(A)                      # = I for direct case

# ---- evaluate REML ll on a grid ---------------------------------------------
rho_grid     <- seq(0.01, 0.99, length.out = 40)
log_phi_grid <- seq(log(0.1), log(1000), length.out = 40)

cat("Evaluating REML surface on 40x40 grid...\n")

ll_mat <- matrix(NA_real_, length(rho_grid), length(log_phi_grid))

for (i in seq_along(rho_grid)) {
  rho_i  <- rho_grid[i]
  prior  <- Q_fun(c(rho = rho_i))
  Q_mat  <- prior$Q_matrix
  ld_Q   <- prior$log_det_Q
  apply_Q      <- function(v) as.numeric(Q_mat %*% v)
  apply_AtRinvA <- function(v) as.numeric(AtRinvA %*% v)

  for (j in seq_along(log_phi_grid)) {
    phi_j    <- exp(log_phi_grid[j])
    apply_K  <- function(v) apply_AtRinvA(v) + (1/phi_j) * apply_Q(v)

    pcg_res  <- tryCatch(
      pcg(apply_K, AtRinvy, tol = 1e-6, maxit = 4L * p),
      error = function(e) list(converged = FALSE)
    )
    if (!pcg_res$converged) next

    yHinvy <- yRinvy - as.numeric(crossprod(AtRinvy, pcg_res$x))
    if (yHinvy <= 0) next
    sigma2e_hat <- yHinvy / n

    K_mat   <- Matrix::forceSymmetric(AtRinvA + (1/phi_j) * Q_mat)
    ld_K    <- tryCatch(
      as.numeric(Matrix::determinant(K_mat, logarithm = TRUE)$modulus),
      error = function(e) NA_real_
    )
    if (!is.finite(ld_K)) next

    ll_mat[i, j] <- -n/2 * log(sigma2e_hat) -
      1/2 * ld_K             -
      p/2 * log(phi_j)       +
      1/2 * ld_Q
  }
}

cat("Done.\n\n")

# ---- find grid maximum -------------------------------------------------------
idx_max  <- which(ll_mat == max(ll_mat, na.rm = TRUE), arr.ind = TRUE)
rho_max  <- rho_grid[idx_max[1]]
phi_max  <- exp(log_phi_grid[idx_max[2]])
cat(sprintf("Grid maximum:  rho = %.3f,  phi = %.3f,  ll = %.3f\n",
            rho_max, phi_max, max(ll_mat, na.rm = TRUE)))
cat(sprintf("True values:   rho = %.3f,  phi = %.3f\n\n", true_rho, phi))

# ---- plot surface ------------------------------------------------------------
# Clip extreme values for cleaner visualisation
ll_plot <- ll_mat
ll_plot[ll_plot < quantile(ll_mat, 0.05, na.rm = TRUE)] <- NA

filled.contour(
  x    = rho_grid,
  y    = log_phi_grid,
  z    = ll_plot,
  xlab = "rho",
  ylab = "log(phi)",
  main = "REML log-likelihood surface\n(A=I, true rho=0.3, true phi=5)",
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
  plot.axes = {
    axis(1); axis(2)
    # mark true values
    abline(v = true_rho,    col = "blue",  lty = 2, lwd = 2)
    abline(h = log(phi),    col = "blue",  lty = 2, lwd = 2)
    # mark grid optimum
    abline(v = rho_max,     col = "green", lty = 1, lwd = 2)
    abline(h = log(phi_max),col = "green", lty = 1, lwd = 2)
    legend("topleft",
           legend = c("true", "grid optimum"),
           col    = c("blue", "green"),
           lty    = c(2, 1), lwd = 2, bg = "white")
  }
)

# ---- also plot 1D slices -----------------------------------------------------
par(mfrow = c(1, 2))

# Slice over rho at true log(phi)
j_true <- which.min(abs(log_phi_grid - log(phi)))
plot(rho_grid, ll_mat[, j_true],
     type = "l", lwd = 2,
     xlab = "rho", ylab = "REML ll",
     main = sprintf("ll vs rho  |  phi fixed at %.1f (true)", phi))
abline(v = true_rho, col = "blue", lty = 2, lwd = 2)

# Slice over log(phi) at true rho
i_true <- which.min(abs(rho_grid - true_rho))
plot(exp(log_phi_grid), ll_mat[i_true, ],
     type = "l", lwd = 2, log = "x",
     xlab = "phi (log scale)", ylab = "REML ll",
     main = sprintf("ll vs phi  |  rho fixed at %.1f (true)", true_rho))
abline(v = phi, col = "blue", lty = 2, lwd = 2)

par(mfrow = c(1, 1))
