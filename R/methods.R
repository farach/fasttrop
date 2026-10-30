#' Methods for `trop` objects
#'
#' Print, summarize, and extract results from a fitted [trop()] model.
#'
#' @param x,object A `trop` object.
#' @param digits Number of significant digits to print.
#' @param parm Ignored; the model has one parameter (the ATT).
#' @param level Confidence level. Only the level used at fit time is
#'   available, because the interval is stored with the fit.
#' @param ... Unused.
#' @name trop-methods
#' @return `coef()` returns the named ATT; `confint()` a 1 x 2 matrix;
#'   `vcov()` a 1 x 1 matrix; `as.data.frame()` one row per treated cell with
#'   its estimated effect; `print()` and `summary()` return `x` or a summary
#'   object invisibly or for printing.
NULL

#' @rdname trop-methods
#' @export
print.trop <- function(x, digits = 4L, ...) {
  cat("Triply Robust Panel (TROP) estimate\n")
  cat(sprintf("  Method: %s | Design: %s | Units: %d | Periods: %d | Treated cells: %d\n",
              x$method, x$design, x$n_units, x$n_periods, x$n_treated_cells))
  cat(sprintf("  ATT: %s", format(signif(x$estimate, digits))))
  if (is.finite(x$std.error)) {
    cat(sprintf("  (SE %s; %g%% CI %s to %s)", format(signif(x$std.error, digits)),
                100 * x$conf.level, format(signif(x$conf.low, digits)),
                format(signif(x$conf.high, digits))))
  }
  cat("\n")
  cat(sprintf("  lambda: time = %s, unit = %s, nn = %s\n",
              format(x$lambda[["time"]]), format(x$lambda[["unit"]]), format(x$lambda[["nn"]])))
  invisible(x)
}

#' @rdname trop-methods
#' @export
summary.trop <- function(object, ...) {
  structure(object, class = c("summary.trop", "trop"))
}

#' @rdname trop-methods
#' @export
print.summary.trop <- function(x, digits = 4L, ...) {
  print.trop(x, digits = digits)
  if (!is.null(x$cv)) {
    cat(sprintf("  Tuning: %s, score %s, %d grid evaluations\n", x$cv$criterion,
                format(signif(x$cv$score, digits)), x$cv$n_evals))
  } else {
    cat("  Tuning: fixed by the user\n")
  }
  if (!is.null(x$boot)) {
    cat(sprintf("  Bootstrap: %d draws (%s), %d failed, %s interval\n", length(x$boot$draws),
                x$boot$resample, x$boot$n_failed, x$ci_type))
  }
  cat(sprintf("  Treated units: %d | Skipped cells: %d | Non-converged fits: %d\n",
              x$n_treated_units, x$n_skipped, x$n_nonconverged))
  cat(sprintf("  Time: %.3f s total, %.3f s in the native kernel (%d threads)\n",
              x$timing$total_seconds, x$timing$kernel_seconds, x$control$threads))
  if (x$design != "block") {
    cat("  Note: Theorem 5.1 (triple robustness) is proven for block assignment only.\n")
  }
  invisible(x)
}

#' @rdname trop-methods
#' @export
coef.trop <- function(object, ...) c(ATT = object$estimate)

#' @rdname trop-methods
#' @export
confint.trop <- function(object, parm, level = object$conf.level, ...) {
  if (!isTRUE(all.equal(level, object$conf.level))) {
    stop("Only the bootstrap interval at the fitted `conf_level` is available; refit to change it.",
         call. = FALSE)
  }
  a <- (1 - level) / 2
  matrix(c(object$conf.low, object$conf.high), 1L, 2L,
         dimnames = list("ATT", sprintf("%g %%", 100 * c(a, 1 - a))))
}

#' @rdname trop-methods
#' @export
vcov.trop <- function(object, ...) matrix(object$std.error^2, 1L, 1L, dimnames = list("ATT", "ATT"))

#' @rdname trop-methods
#' @param row.names,optional Ignored.
#' @export
as.data.frame.trop <- function(x, row.names = NULL, optional = FALSE, ...) {
  idx <- which(!is.na(x$tau), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  units <- if (!is.null(x$units)) x$units[idx[, 1]] else idx[, 1]
  times <- if (!is.null(x$times)) x$times[idx[, 2]] else idx[, 2]
  data.frame(unit = units, time = times, tau = x$tau[idx], row.names = NULL)
}
