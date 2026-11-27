#' Count the cells of each group before the bootstrap
#'
#' `trop_cells()` shows which treated cells are in each group of
#' [trop_effects()], and which cells are not. It does not show the
#' effects. Thus, you can examine the estimand before you see the results or
#' calculate the bootstrap.
#'
#' The function uses the groups, the window, the anticipation, and the weights
#' of [trop_effects()]. For each group, it counts the units and the cells in
#' the estimate. It also counts the treated cells outside of the target of the
#' fit, and the cells that the untreated cells do not identify. A group with
#' only cells outside of the estimate is also in the result.
#'
#' The column `share` is the weight of each group in the overall estimate:
#' the mean of all cells in the window with the same weights. With equal
#' weights for each cell, the share of a group is its number of cells divided
#' by the number of cells in the estimate.
#'
#' If no treated cell is in a group, the result has no rows. The attribute
#' `totals` and the print method then show where the treated cells are.
#'
#' @inheritParams trop_effects
#'
#' @return A data frame of class `trop_cells` with one row for each group and
#'   the columns `n_units`, `n_cells`, `n_excluded`, `n_skipped`, and `share`.
#'   The column `n_excluded` counts the treated cells with an outcome that are
#'   outside of the target. The column `n_skipped` counts the cells in the
#'   target that have no effect, because the untreated cells do not identify
#'   them or their fit failed. The attribute `totals` has the counts for all
#'   treated cells.
#' @seealso [trop_effects()], [trop_rollout()]
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 6,
#'                      n_post = 4, design = "staggered", seed = 1)
#' sim$in_target <- sim$time <= 11
#' fit <- suppressWarnings(
#'   trop(sim, "y", "d", "unit", "time", lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'        target = "in_target", se = "none", control = trop_control(threads = 1))
#' )
#' trop_cells(fit, by = "event_time")
#' @export
trop_cells <- function(fit, by = c("event_time", "cohort", "time", "unit", "overall"),
                       window = NULL, anticipation = 0L, weights = c("cell", "unit"),
                       unit_weights = NULL) {
  check_local_fit(fit, "trop_cells")
  weights <- match.arg(weights)
  anticipation <- check_anticipation(anticipation)
  check_event_source(fit, anticipation)
  if (is.list(by) && !is.data.frame(by)) {
    stop("`trop_cells()` takes one set of groups in `by`, not a list.", call. = FALSE)
  }
  P <- fit$panel
  W <- P$W
  N <- nrow(W)
  TT <- ncol(W)
  wt <- effect_weights(weights, unit_weights, fit)
  excluded <- if (is.null(P$excluded)) matrix(FALSE, N, TT) else P$excluded
  has_y <- is.finite(P$Y)
  treated <- W == 1L & (has_y | excluded)
  in_est <- W == 1L & is.finite(fit$tau) & wt$keep[row(W)]
  # All treated cells with an outcome get a group, so the cells outside of the
  # estimate are also in the groups.
  mark <- matrix(NA_real_, N, TT)
  mark[treated] <- 0
  if (!is.null(window) || (is.character(by) && by[1L] %in% c("event_time", "cohort"))) {
    anticipation_warning(W, anticipation, P$event)
  }
  grp <- effect_groups(by, W, mark, window, fit$units, fit$times, anticipation,
                       event_given = P$event, allow_empty = TRUE)
  id <- grp$id
  G <- length(grp$labels)
  in_g <- !is.na(id)
  est <- in_g & in_est
  skipped <- in_g & W == 1L & has_y & !is.finite(fit$tau)
  ex <- in_g & excluded
  check_weight_range(wt$u, list(ifelse(est, id, NA_integer_)))
  ui <- row(W)[est]
  gi <- id[est]
  u <- if (is.null(wt$u)) rep(1, length(ui)) else wt$u[ui]
  # The shares compare the groups. Thus, one scale for all cells.
  if (!is.null(wt$u) && length(u) > 0L) u <- u / max(u)
  if (wt$per_unit && length(ui) > 0L) u <- u / stats::ave(ui, ui, FUN = length)
  share <- numeric(G)
  if (length(gi) > 0L) {
    s <- rowsum(u, gi)
    share[as.integer(rownames(s))] <- s[, 1L] / sum(u)
  }
  pairs <- !duplicated(cbind(ui, gi))
  out <- data.frame(group = grp$labels, n_units = tabulate(gi[pairs], G),
                    n_cells = tabulate(gi, G), n_excluded = tabulate(id[ex], G),
                    n_skipped = tabulate(id[skipped], G), share = share)
  names(out)[1L] <- grp$name
  all_grp <- if (is.null(window)) grp else {
    effect_groups(by, W, mark, NULL, fit$units, fit$times, anticipation,
                  event_given = P$event, allow_empty = TRUE)
  }
  in_any <- !is.na(all_grp$id)
  # The treated cells with an outcome are in exactly one of these counts:
  # in_estimate, excluded, skipped, zero_weight, outside_window, and no_group.
  attr(out, "totals") <- c(
    treated = sum(W == 1L), with_outcome = sum(treated), in_estimate = sum(est),
    excluded = sum(ex), skipped = sum(skipped),
    zero_weight = sum(in_g & W == 1L & is.finite(fit$tau) & !wt$keep[row(W)]),
    outside_window = sum(in_any) - sum(in_g), no_group = sum(treated) - sum(in_any),
    no_outcome = sum(W == 1L & !treated)
  )
  attr(out, "by") <- grp$name
  attr(out, "anticipation") <- anticipation
  attr(out, "window") <- window
  attr(out, "weights") <- wt$label
  attr(out, "event_given") <- !is.null(P$event)
  class(out) <- c("trop_cells", "data.frame")
  out
}

#' @rdname trop_cells
#' @param x A `trop_cells` object.
#' @param digits The number of significant digits to show.
#' @param ... Not used.
#' @export
print.trop_cells <- function(x, digits = 3L, ...) {
  cat(sprintf("TROP cells %s\n", group_title(attr(x, "by"))))
  for (line in effects_notes(x)) cat("  ", line, "\n", sep = "")
  df <- x
  class(df) <- "data.frame"
  for (a in c("totals", "by", "anticipation", "window", "weights", "event_given")) attr(df, a) <- NULL
  if (nrow(df) > 0L) {
    print(df, digits = digits, row.names = FALSE)
  } else {
    cat("No treated cell is in a group.\n")
  }
  tot <- attr(x, "totals")
  cat(sprintf("Cells in the estimate: %d. The column share is the weight of each group in the overall estimate.\n",
              tot[["in_estimate"]]))
  other <- c(
    if (tot[["excluded"]] > 0L) sprintf("%d outside of the target", tot[["excluded"]]),
    if (tot[["skipped"]] > 0L) sprintf("%d with no effect", tot[["skipped"]]),
    if (tot[["zero_weight"]] > 0L) sprintf("%d with a weight of 0", tot[["zero_weight"]]),
    if (tot[["outside_window"]] > 0L) sprintf("%d outside of the window", tot[["outside_window"]]),
    if (tot[["no_group"]] > 0L) sprintf("%d in no group", tot[["no_group"]]),
    if (tot[["no_outcome"]] > 0L) sprintf("%d with no outcome", tot[["no_outcome"]])
  )
  if (length(other) > 0L) {
    cat(sprintf("Treated cells that are not in the estimate: %s.\n", paste(other, collapse = ", ")))
  }
  invisible(x)
}
