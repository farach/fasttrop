#' Weighted fit for given weights
#'
#' Internal helper. Minimizes
#' `sum(w * (Y - alpha_i - beta_t - L_it)^2) + lambda_nn * ||L||_*`
#' with the package's native solver. `lambda_nn = Inf` drops the low-rank term.
#' Used by the tests to compare the solver with reference implementations on
#' identical weights.
#' @noRd
trop_fit_weighted <- function(Y, w, lambda_nn, control = trop_control()) {
  if (!is.matrix(Y) || !is.matrix(w) || !all(dim(Y) == dim(w))) {
    stop("`Y` and `w` must be matrices with the same dimensions.", call. = FALSE)
  }
  if (anyNA(w) || any(w < 0)) stop("`w` must be non-negative without missing values.", call. = FALSE)
  if (!inherits(control, "trop_control")) control <- do.call(trop_control, as.list(control))
  Y[w == 0] <- 0
  if (anyNA(Y)) stop("`Y` must be finite wherever `w` is positive.", call. = FALSE)
  storage.mode(Y) <- "double"
  storage.mode(w) <- "double"
  .Call(C_fit_weighted, Y, w, as.double(lambda_nn), unclass(control))
}
