#' Triply Robust Panel (TROP) estimator
#'
#' Estimates the average treatment effect on the treated (ATT) with the TROP
#' estimator of Athey, Imbens, Qu, and Viviano (2026). For each treated cell,
#' the counterfactual untreated outcome comes from a weighted two-way fixed
#' effects model with a nuclear-norm penalized low-rank term, fitted on
#' untreated cells. Unit weights decay with the distance between pre-treatment
#' outcome paths; time weights decay with the distance in periods.
#'
#' @section Methods:
#' * `method = "local"` follows the paper's Eq. 12 and Algorithm 2. Every
#'   treated cell gets its own unit and time weights and its own fit. Tuning
#'   uses the leave-one-out criterion of Eq. 4-5 over untreated cells and the
#'   coordinate search of the paper's footnote 2. Works for block, staggered,
#'   and general assignment.
#' * `method = "global"` matches the authors' Python reference package: one
#'   fit for a block design, with time distance to the centre of the treated
#'   block and unit distance to the average treated pre-period path. Tuning uses
#'   placebo cross-validation on untreated units (coordinate cycling from the
#'   grid means). Requires a block design.
#'
#' @section Inference:
#' Standard errors come from the unit bootstrap of the paper's Algorithm 3,
#' holding the tuning parameters fixed. `resample = "stratified"` resamples
#' treated and untreated units separately; `"controls"` keeps the treated units
#' fixed. The standard error is the standard deviation of the draws. The default
#' interval is the normal interval built from it (`ci = "normal"`), which
#' follows the paper's bootstrap variance. `ci = "percentile"` uses quantiles of
#' the draws, as the authors' Python package does; it needs more draws for
#' stable tail quantiles. The paper notes that bootstrap validity requires a
#' growing number of treated units, so the package warns when fewer than 5
#' units are treated.
#'
#' @section Scope of the triple-robustness result:
#' The paper's Theorem 5.1 (bias bounded by the product of unit imbalance, time
#' imbalance, and regression-adjustment error) is proven under block
#' assignment. For staggered or general designs the estimator is defined, but
#' that guarantee is not established, and the package warns.
#'
#' @param data A data frame in long format, one row per unit-period.
#' @param outcome,treatment,unit,time Column names (strings). Treatment must be
#'   0/1. Missing outcomes are allowed and are excluded from the fits.
#' @param method `"local"` or `"global"`; see Methods.
#' @param lambda Optional tuning parameters as a named numeric vector or list
#'   with `time`, `unit`, and `nn`. If `NULL`, they are chosen by
#'   cross-validation over `grid`.
#' @param grid Tuning grids from [trop_grid()].
#' @param cv_cells For `method = "local"`: number of untreated cells sampled
#'   for the leave-one-out criterion. `NULL` uses every untreated cell (Eq. 5).
#' @param n_placebo,placebo_size For `method = "global"`: number of placebo
#'   sets and untreated units per set (defaults follow the reference package:
#'   200 sets of 1 unit).
#' @param max_cycles Maximum coordinate-search cycles.
#' @param se `"bootstrap"` or `"none"`.
#' @param n_boot Number of bootstrap draws.
#' @param resample `"stratified"` or `"controls"`; see Inference.
#' @param conf_level Confidence level of the interval.
#' @param ci `"normal"` (estimate plus or minus a normal quantile times the
#'   bootstrap standard error) or `"percentile"` (quantiles of the bootstrap
#'   draws); see Inference.
#' @param control Solver settings from [trop_control()].
#' @param seed Optional integer seed for cell sampling, placebo sets, and the
#'   bootstrap. The global random state is restored afterwards.
#'
#' @return An object of class `trop` with the ATT (`estimate`), bootstrap
#'   standard error and confidence interval, cell-level effects (`tau`), the
#'   tuning parameters, the cross-validation trace, convergence diagnostics,
#'   and timings that separate native kernel time from total time.
#'
#' @references Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026).
#'   Triply Robust Panel Estimators. *Journal of Applied Econometrics*.
#'   \doi{10.1002/jae.70061}.
#'
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' fit <- trop(sim, "y", "d", "unit", "time",
#'             lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'             n_boot = 50, seed = 1, control = trop_control(threads = 1))
#' fit
#' @export
trop <- function(data, outcome, treatment, unit, time,
                 method = c("local", "global"), lambda = NULL, grid = trop_grid(),
                 cv_cells = NULL, n_placebo = 200L, placebo_size = 1L, max_cycles = 10L,
                 se = c("bootstrap", "none"), n_boot = 200L,
                 resample = c("stratified", "controls"), conf_level = 0.95,
                 ci = c("normal", "percentile"), control = trop_control(), seed = NULL) {
  t0 <- Sys.time()
  pm <- panel_matrices(data, outcome, treatment, unit, time)
  fit <- trop_matrix(pm$Y, pm$W, method = method, lambda = lambda, grid = grid,
                     cv_cells = cv_cells, n_placebo = n_placebo, placebo_size = placebo_size,
                     max_cycles = max_cycles, se = se, n_boot = n_boot, resample = resample,
                     conf_level = conf_level, ci = ci, control = control, seed = seed)
  fit$units <- pm$units
  fit$times <- pm$times
  fit$call <- match.call()
  fit$timing$total_seconds <- as.double(difftime(Sys.time(), t0, units = "secs"))
  fit
}
