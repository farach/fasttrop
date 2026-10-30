#' Weighted fit for given weights
#'
#' This function is internal. It finds the minimum of
#' `sum(w * (Y - alpha_i - beta_t - L_it)^2) + lambda_nn * ||L||_*`
#' with the C solver. `lambda_nn = Inf` removes the low-rank term. The tests
#' use this function to compare the solver with the reference implementations
#' on the same weights.
#' @noRd
trop_fit_weighted <- function(Y, w, lambda_nn, control = trop_control()) {
  if (!is.matrix(Y) || !is.matrix(w) || !all(dim(Y) == dim(w))) {
    stop("`Y` and `w` must be matrices with the same dimensions.", call. = FALSE)
  }
  if (anyNA(w) || any(w < 0)) {
    stop("`w` must contain only values of 0 or more and no `NA` values.", call. = FALSE)
  }
  if (!inherits(control, "trop_control")) control <- do.call(trop_control, as.list(control))
  Y[w == 0] <- 0
  if (anyNA(Y)) stop("`Y` must be finite in each cell where `w` is more than 0.", call. = FALSE)
  storage.mode(Y) <- "double"
  storage.mode(w) <- "double"
  .Call(C_fit_weighted, Y, w, as.double(lambda_nn), unclass(control))
}
