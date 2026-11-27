# Diagnostics of the fits of the bootstrap and of the cross-validation.

# A data frame with the counts of the fits of each draw, from the C code. The
# column n_dropped counts the cells of the draw that have an effect in the fit
# on the data, but that the untreated cells of the draw do not identify.
boot_diagnostics <- function(m) {
  m <- matrix(as.integer(m), ncol = 5L)
  data.frame(draw = seq_len(nrow(m)), n_used = m[, 1L], n_skipped = m[, 2L],
             n_nonconverged = m[, 3L], n_failed = m[, 4L], n_dropped = m[, 5L])
}

boot_warnings <- function(d) {
  if (is.null(d)) return(invisible(NULL))
  B <- nrow(d)
  if (any(d$n_failed > 0L)) {
    warning(sprintf("In %d of %d bootstrap draws, %d fit(s) failed: the singular value decomposition failed or the result is not finite. These draws do not use the cells of these fits.",
                    sum(d$n_failed > 0L), B, sum(d$n_failed)), call. = FALSE)
  }
  if (any(d$n_nonconverged > 0L)) {
    warning(sprintf("In %d of %d bootstrap draws, %d fit(s) did not converge. Use a larger `max_iter` in trop_control().",
                    sum(d$n_nonconverged > 0L), B, sum(d$n_nonconverged)), call. = FALSE)
  }
  if (any(d$n_dropped > 0L)) {
    warning(sprintf("In %d of %d bootstrap draws, the untreated cells do not identify %d treated cell(s) that have an effect in the fit. These draws use fewer cells.",
                    sum(d$n_dropped > 0L), B, sum(d$n_dropped)), call. = FALSE)
  }
  invisible(NULL)
}

# Short lines for the print methods. A line for each kind of problem.
diag_lines <- function(d) {
  if (is.null(d)) return(character())
  out <- character()
  if (any(d$n_nonconverged > 0L)) {
    out <- c(out, sprintf("%d fit(s) in %d draw(s) did not converge.",
                          sum(d$n_nonconverged), sum(d$n_nonconverged > 0L)))
  }
  if (any(d$n_failed > 0L)) {
    out <- c(out, sprintf("%d fit(s) in %d draw(s) failed.", sum(d$n_failed), sum(d$n_failed > 0L)))
  }
  if (any(d$n_dropped > 0L)) {
    out <- c(out, sprintf("%d treated cell(s) in %d draw(s) have an effect in the fit, but not in the draw.",
                          sum(d$n_dropped), sum(d$n_dropped > 0L)))
  }
  out
}

# Format a time in seconds for a message.
format_seconds <- function(s) {
  if (!is.finite(s)) return("unknown")
  if (s < 90) return(sprintf("%.0f s", s))
  if (s < 5400) return(sprintf("%.1f min", s / 60))
  sprintf("%.1f h", s / 3600)
}
