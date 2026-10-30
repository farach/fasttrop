#' fasttrop: Fast Triply Robust Panel Estimators
#'
#' The package calculates the Triply Robust Panel (TROP) estimator of Athey,
#' Imbens, Qu, and Viviano (2026). The estimator adds unit weights, time
#' weights, and a low-rank term with a nuclear-norm penalty to a model with
#' two-way fixed effects. All calculations run in C. The C code can use more
#' than one thread with 'OpenMP'.
#'
#' Use [trop()] for a data frame in long format. Use [trop_matrix()] for a
#' matrix of outcomes and a matrix of treatment indicators.
#'
#' @references Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026).
#'   Triply Robust Panel Estimators. *Journal of Applied Econometrics*.
#'   \doi{10.1002/jae.70061}. Preprint: arXiv:2508.21536.
#'
#' @keywords internal
#' @useDynLib fasttrop, .registration = TRUE
"_PACKAGE"
