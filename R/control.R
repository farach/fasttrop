#' Solver options
#'
#' `trop_control()` sets the options of the weighted low-rank solver in C. The
#' estimate, the cross-validation, and the bootstrap all use this solver.
#'
#' @param max_iter The maximum number of outer iterations. With
#'   `solver = "block"`, each outer iteration is one sweep of the blocks. With
#'   `solver = "prox"`, the maximum number of proximal steps is
#'   `max_iter * max_inner`.
#' @param max_inner The number of proximal steps in each outer iteration. Each
#'   step is a step of FISTA, an accelerated proximal gradient method.
#' @param tol The convergence tolerance. The solver stops when the largest
#'   absolute change of the iterates between two steps is less than `tol`.
#' @param threads The number of 'OpenMP' threads. The default comes from
#'   [trop_threads()].
#' @param warm Logical. If `TRUE`, the fit for a cell starts from the solution
#'   for the previous cell of the same unit. If `FALSE`, each fit starts from
#'   zero. Both options give the same results in the limits of `tol`. The
#'   number of threads does not change the results.
#' @param solver Use `"prox"`, the default, or `"block"`. `"prox"` is an
#'   accelerated proximal gradient method on the low-rank term. At each step, it
#'   calculates the unit and time effects again, and it uses an adaptive
#'   restart. `"block"` is the block coordinate method of the reference
#'   implementations. Both methods find the minimum of the same convex
#'   objective. `"prox"` usually calculates fewer singular value decompositions.
#' @param throttle Logical. Windows can use power throttling (EcoQoS) to make
#'   processes that run in the background for a long time slower. An R session
#'   with a long calculation can be such a process.
#'
#'   With `FALSE`, the default, each thread that runs the C code of fasttrop
#'   tells Windows not to use power throttling on it. When the C code stops,
#'   Windows controls the thread again. With `TRUE`, Windows uses its usual
#'   policy. On other operating systems, this option has no effect.
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
    stop("`tol` must be one positive number.", call. = FALSE)
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
#' `trop_threads()` gives the default number of threads for [trop_control()].
#'
#' @return An integer. If the option `fasttrop.threads` has a value, the
#'   function returns this value. If not, the function returns the number of
#'   physical cores, with a minimum of 1. When the environment variable
#'   `_R_CHECK_LIMIT_CORES_` has a value other than `"false"`, the maximum is 2.
#'   `R CMD check --as-cran` sets this variable.
#' @examples
#' trop_threads()
#' @export
trop_threads <- function() {
  opt <- getOption("fasttrop.threads")
  if (!is.null(opt)) return(as.integer(opt))
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(n) || n < 1L) n <- 1L
  # CRAN policy: package checks must not use more than 2 cores.
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && limit != "false") n <- min(n, 2L)
  as.integer(n)
}

#' Tuning grids
#'
#' `trop_grid()` gives the candidate values of the three tuning parameters.
#' These parameters are the decay of the time weights, the decay of the unit
#' weights, and the nuclear-norm penalty. A value of 0 for `time` or `unit` gives
#' equal weights. A value of `Inf` for `nn` removes the low-rank term. Then the
#' model has only two-way fixed effects.
#'
#' @param time The candidate values for the decay of the time weights. The
#'   values must be finite and 0 or more.
#' @param unit The candidate values for the decay of the unit weights. The
#'   values must be finite and 0 or more.
#' @param nn The candidate values for the nuclear-norm penalty. The values must
#'   be 0 or more. You can use `Inf`.
#'
#' @return A list of class `trop_grid` with the numeric vectors `time`, `unit`,
#'   and `nn`.
#' @examples
#' trop_grid(nn = c(0.1, 1, Inf))
#' @export
trop_grid <- function(time = c(0, 0.1, 0.5, 1, 2, 5),
                      unit = c(0, 0.1, 0.5, 1, 2, 5),
                      nn = c(0.01, 0.1, 1, 10, Inf)) {
  ok <- function(x, finite) {
    is.numeric(x) && length(x) >= 1L && !anyNA(x) && all(x >= 0) && (!finite || all(is.finite(x)))
  }
  if (!ok(time, TRUE)) stop("`time` must contain only finite numbers of 0 or more.", call. = FALSE)
  if (!ok(unit, TRUE)) stop("`unit` must contain only finite numbers of 0 or more.", call. = FALSE)
  if (!ok(nn, FALSE)) stop("`nn` must contain only numbers of 0 or more. You can use `Inf`.", call. = FALSE)
  structure(list(time = as.double(time), unit = as.double(unit), nn = as.double(nn)),
            class = "trop_grid")
}

check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x)) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
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
    stop("In `lambda`, `time` and `unit` must be finite numbers of 0 or more, and `nn` must be 0 or more.",
         call. = FALSE)
  }
  lambda
}
