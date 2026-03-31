#' Predict from a fastblm_fit object
#'
#' @param object fastblm_fit object
#' @param A_new n_new x p design matrix for prediction
#' @param ... unused
#'
#' @return numeric vector of predicted means
#' @export
predict.fastblm_fit <- function(object, A_new, ...) {
  as.numeric(A_new %*% object$posterior_mean)
}
