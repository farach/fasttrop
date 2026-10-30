#' fasttrop: Fast Triply Robust Panel Estimators
#'
#' Implements the Triply Robust Panel (TROP) estimator of Athey, Imbens, Qu,
#' and Viviano (2026), which combines unit weights, time weights, and a
#' nuclear-norm regularized factor model on top of two-way fixed effects.
#' The numerical work (weighted low-rank fits, cross-validation, and the
#' bootstrap) runs in C with optional 'OpenMP' threads.
#'
#' Main entry points are [trop()] for long-format data frames and
#' [trop_matrix()] for outcome and treatment matrices.
#'
#' @references Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026).
#'   Triply Robust Panel Estimators. *Journal of Applied Econometrics*.
#'   \doi{10.1002/jae.70061}. Preprint: arXiv:2508.21536.
#'
#' @keywords internal
#' @useDynLib fasttrop, .registration = TRUE
"_PACKAGE"
