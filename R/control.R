#' Solver controls
#'
#' Options for the native weighted low-rank solver shared by estimation,
#' cross-validation, and the bootstrap.
#'
#' @param max_iter Maximum outer iterations. For `solver = "block"` these are
#'   block-coordinate sweeps; for `solver = "prox"` the step budget is
#'   `max_iter * max_inner` accelerated proximal steps.
#' @param max_inner Proximal (FISTA) iterations per outer iteration.
#' @param tol Convergence tolerance on the maximum absolute change of the
#'   iterates between steps.
#' @param threads Number of 'OpenMP' threads. Defaults to the option
#'   `fasttrop.threads`, or the number of physical cores.
#' @param warm Logical. If `TRUE`, fits for successive cells of the same unit
#'   start from the previous solution. Results then agree with cold starts up
#'   to the convergence tolerance and remain deterministic for any thread count.
#' @param solver `"prox"` (default) runs an accelerated proximal gradient on the
#'   low-rank term with the unit and time effects profiled out at every step
#'   and adaptive restart. `"block"` runs the block-coordinate scheme of the
#'   reference implementations (one fixed-effect sweep, then `max_inner`
#'   proximal steps per outer iteration). Both minimize the same convex
#'   objective and agree at convergence; `"prox"` usually needs far fewer
#'   singular value decompositions.
#' @param throttle Logical. Windows may slow down long-running background
#'   processes to save power (EcoQoS power throttling). In a measured R session
#'   this made identical native work about 1.3 times slower than in a fresh
#'   process. With the default `FALSE`, each thread running a fasttrop kernel
#'   asks Windows not to throttle it while the kernel runs, then hands the
#'   decision back to the system. Use `TRUE` to keep the system policy.
#'   Ignored on other platforms.
#'
#' @return A list of class `trop_control`.
#' @examples
#' trop_control(tol = 1e-8, threads = 2)
#' @export
trop_control <- function(max_iter = 100L, max_inner = 10L, tol = 1e-6,
                         threads = trop_threads(), warm = TRUE,
                         solver = c("prox", "block"), throttle = FALSE) {
  solver <- match.arg(solver)
  check_count(max_iter, "max_iter")
  check_count(max_inner, "max_inner")
  check_count(threads, "threads")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }
  if (!is.logical(warm) || length(warm) != 1L || is.na(warm)) {
    stop("`warm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(throttle) || length(throttle) != 1L || is.na(throttle)) {
    stop("`throttle` must be TRUE or FALSE.", call. = FALSE)
  }
  structure(
    list(max_iter = as.integer(max_iter), max_inner = as.integer(max_inner),
         tol = as.double(tol), threads = as.integer(threads), warm = as.integer(warm),
         solver = if (solver == "block") 1L else 0L, throttle = as.integer(throttle)),
    class = "trop_control"
  )
}

#' Default number of threads
#'
#' @return The option `fasttrop.threads` if set, otherwise the number of
#'   physical cores (at least 1).
#' @examples
#' trop_threads()
#' @export
trop_threads <- function() {
  opt <- getOption("fasttrop.threads")
  if (!is.null(opt)) return(as.integer(opt))
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(n) || n < 1L) n <- 1L
  as.integer(n)
}

#' Tuning grids
#'
#' Candidate values for the time-decay, unit-decay, and nuclear-norm tuning
#' parameters. A value of 0 for `time` or `unit` gives uniform weights;
#' `Inf` for `nn` disables the low-rank component (two-way fixed effects only).
#'
#' @param time Candidate values of the time-weight decay `lambda_time` (finite, >= 0).
#' @param unit Candidate values of the unit-weight decay `lambda_unit` (finite, >= 0).
#' @param nn Candidate values of the nuclear-norm penalty `lambda_nn` (> 0, may include `Inf`).
#'
#' @return A list with numeric vectors `time`, `unit`, and `nn`.
#' @examples
#' trop_grid(nn = c(0.1, 1, Inf))
#' @export
trop_grid <- function(time = c(0, 0.1, 0.5, 1, 2, 5),
                      unit = c(0, 0.1, 0.5, 1, 2, 5),
                      nn = c(0.01, 0.1, 1, 10, Inf)) {
  ok <- function(x, finite) {
    is.numeric(x) && length(x) >= 1L && !anyNA(x) && all(x >= 0) && (!finite || all(is.finite(x)))
  }
  if (!ok(time, TRUE)) stop("`time` must be finite non-negative numbers.", call. = FALSE)
  if (!ok(unit, TRUE)) stop("`unit` must be finite non-negative numbers.", call. = FALSE)
  if (!ok(nn, FALSE)) stop("`nn` must be non-negative numbers (Inf allowed).", call. = FALSE)
  structure(list(time = as.double(time), unit = as.double(unit), nn = as.double(nn)),
            class = "trop_grid")
}

check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x)) {
    stop(sprintf("`%s` must be a single positive integer.", name), call. = FALSE)
  }
  invisible(x)
}

as_lambda <- function(lambda) {
  if (is.list(lambda)) lambda <- unlist(lambda)
  if (!is.numeric(lambda) || length(lambda) != 3L || anyNA(lambda)) {
    stop("`lambda` must be numeric with elements `time`, `unit`, and `nn`.", call. = FALSE)
  }
  if (!is.null(names(lambda)) && all(c("time", "unit", "nn") %in% names(lambda))) {
    lambda <- lambda[c("time", "unit", "nn")]
  }
  lambda <- stats::setNames(as.double(lambda), c("time", "unit", "nn"))
  if (any(lambda < 0) || !is.finite(lambda[["time"]]) || !is.finite(lambda[["unit"]])) {
    stop("`lambda` must have finite non-negative `time` and `unit` and non-negative `nn`.",
         call. = FALSE)
  }
  lambda
}
