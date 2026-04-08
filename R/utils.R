# NULL coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Error if x is a function -- used in Cholesky path where we need explicit matrices
stop_if_function <- function(x, name = "x") {
  if (is.function(x))
    stop(sprintf("`%s` is a function -- supply an explicit matrix or use solver = 'pcg'.", name))
  invisible(x)
}

# Coerce matrix or function to apply function -- used in PCG path
as_apply <- function(x) {
  if (is.null(x))     return(function(v) v)   # identity
  if (is.function(x)) return(x)
  function(v) as.numeric(x %*% v)
}

# Resolve R_inv: NULL means identity
resolve_Rinv <- function(R_inv, n) {
  if (is.null(R_inv)) return(Matrix::Diagonal(n))
  R_inv
}

# Check that an object is a matrix (base or Matrix)
is_matrix <- function(x) is.matrix(x) || inherits(x, "Matrix")

# Assemble the posterior precision operator K(phi)
# K(phi) v = A' R_inv A v + (1/phi) Q v
make_apply_K <- function(apply_A, apply_At, apply_Q, apply_Rinv, phi) {
  function(v) apply_At(apply_Rinv(apply_A(v))) + (1/phi) * apply_Q(v)
}
