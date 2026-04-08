# =============================================================================

# diagnose_phi_slice.R

#

# For the A=I, rho=0.3 dataset, plots how each term of the REML ll

# varies with phi at the true rho. This reveals whether the flat phi

# surface is a cancellation between terms, or something unexpected.

# =============================================================================

set.seed(42)

m <- 20; p <- m * m; phi_true <- 5.0; sigma2e_true <- 1.0

W        <- make_grid_W(m)

true_rho <- 0.3

x        <- simulate_sar(W, true_rho, phi_true, sigma2e_true)

y        <- x + rnorm(p, sd = sqrt(sigma2e_true))

A        <- diag(p)

Q_fun    <- make_sar_Q_fun(W)

prior    <- Q_fun(c(rho = true_rho))

Q_mat    <- prior$Q_matrix

ld_Q     <- prior$log_det_Q

AtA      <- crossprod(A)           # = I for direct case

AtRinvy  <- as.numeric(crossprod(A, y))

yRinvy   <- as.numeric(crossprod(y))

n        <- length(y)

log_phi_grid <- seq(log(0.1), log(1000), length.out = 200)

# storage

out <- data.frame(

  log_phi   = log_phi_grid,

  phi       = exp(log_phi_grid),

  sigma2e   = NA_real_,

  ll        = NA_real_,

  term_res  = NA_real_,   # -n/2 * log(sigma2e)

  term_logK = NA_real_,   # -1/2 * logdet_K

  term_phi  = NA_real_,   # -p/2 * log(phi)

  term_Q    = NA_real_    # +1/2 * logdet_Q  (constant in phi)

)

for (i in seq_len(nrow(out))) {

  phi_i   <- out$phi[i]

  K_mat   <- Matrix::forceSymmetric(AtA + (1/phi_i) * Q_mat)

  ld_K    <- as.numeric(Matrix::determinant(K_mat, logarithm = TRUE)$modulus)

  x_hat   <- as.numeric(solve(K_mat, AtRinvy))

  yHinvy  <- yRinvy - as.numeric(crossprod(AtRinvy, x_hat))

  if (yHinvy <= 0) next

  s2e     <- yHinvy / n

  out$sigma2e[i]   <- s2e

  out$term_res[i]  <- -n/2 * log(s2e)

  out$term_logK[i] <- -1/2 * ld_K

  out$term_phi[i]  <- -p/2 * log(phi_i)

  out$term_Q[i]    <- 1/2  * ld_Q

  out$ll[i]        <- out$term_res[i] + out$term_logK[i] +

    out$term_phi[i] + out$term_Q[i]

}

# ---- plot -------------------------------------------------------------------

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

plot(out$phi, out$sigma2e, type = "l", lwd = 2, log = "x",

     xlab = "phi", ylab = "sigma2e_hat",

     main = "Profiled sigma2e vs phi")

abline(v = phi_true, col = "blue", lty = 2)

abline(h = sigma2e_true, col = "red", lty = 2)

legend("topright", c("true phi", "true sigma2e"),

       col = c("blue","red"), lty = 2, cex = 0.8)

plot(out$phi, out$term_res, type = "l", lwd = 2, col = "darkgreen", log = "x",

     xlab = "phi", ylab = "value", main = "-n/2 * log(sigma2e)  [residual term]")

abline(v = phi_true, col = "blue", lty = 2)

plot(out$phi, out$term_logK, type = "l", lwd = 2, col = "orange", log = "x",

     xlab = "phi", ylab = "value", main = "-1/2 * log|K|  [complexity term]")

abline(v = phi_true, col = "blue", lty = 2)

plot(out$phi, out$term_phi, type = "l", lwd = 2, col = "purple", log = "x",

     xlab = "phi", ylab = "value", main = "-p/2 * log(phi)  [phi penalty]")

abline(v = phi_true, col = "blue", lty = 2)

# residual + logK combined (the two that move with phi)

plot(out$phi, out$term_res + out$term_logK, type = "l", lwd = 2,

     col = "red", log = "x",

     xlab = "phi", ylab = "value",

     main = "residual + complexity  (sum of moving terms)")

abline(v = phi_true, col = "blue", lty = 2)

plot(out$phi, out$ll, type = "l", lwd = 2, log = "x",

     xlab = "phi", ylab = "REML ll", main = "Total REML ll")

abline(v = phi_true, col = "blue", lty = 2)

par(mfrow = c(1, 1))

# print where optimum is

cat(sprintf("True phi     = %.2f\n", phi_true))

cat(sprintf("Grid optimum = %.2f  (ll = %.4f)\n",

            out$phi[which.max(out$ll)], max(out$ll, na.rm = TRUE)))

cat(sprintf("ll at phi=1000: %.4f\n", out$ll[nrow(out)]))

cat(sprintf("ll drop from optimum to phi=1000: %.4f\n",

            max(out$ll, na.rm=TRUE) - out$ll[nrow(out)]))
