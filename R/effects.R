#' Effects by event time, cohort, period, or unit
#'
#' `trop_effects()` averages the effects of the treated cells of a [trop()] fit
#' in groups. A group is, for example, all treated cells with the same number
#' of periods since the start of the treatment. The function can also give
#' bootstrap standard errors, intervals, and simultaneous bands for all groups.
#'
#' The function uses the effects `tau` of the fit. Thus, the target cells of the
#' fit are the only cells in the groups. A cell that has no effect in the fit,
#' for example a cell outside of the target, is not in a group. A group with no
#' such cell is not in the result. Use [trop_cells()] to count the cells of
#' each group before you calculate the bootstrap.
#'
#' @section Groups:
#' The argument `by` sets the groups:
#' * `"event_time"`: the number of periods since the start of the current
#'   treatment. The value 0 is the first treated period. If a unit stops the
#'   treatment and starts it again, the count starts again.
#' * `"cohort"`: the period in which the current treatment started.
#' * `"time"`: the period.
#' * `"unit"`: the unit.
#' * `"overall"`: one group with all cells. The result is the ATT.
#' * A matrix with the dimensions of the panel. The values are the group of each
#'   cell. A cell with `NA` is not in a group. Only treated cells with an effect
#'   are in groups.
#' * A named list of the values above, for example
#'   `list(event = "event_time", cohort = "cohort")`. The function then gives
#'   a named list with one result for each element. All results use the same
#'   bootstrap draws, and the function calculates the estimator only one time
#'   for each draw. A character element without a name uses its value as the
#'   name.
#'
#' The argument `window` limits the cells to cells with an event time in the
#' window. For example, `window = c(0, 3)` uses the first 4 treated periods.
#'
#' @section Anticipation:
#' Sometimes the treatment has an effect before it starts, for example when
#' units know the start date. A usual procedure is to set `W = 1` in the
#' periods before the start, so that these cells are not controls. Then, the
#' first `W = 1` cell is not the true start.
#'
#' With `anticipation = k`, the
#' function counts the first `k` treated periods of each run as the periods
#' before the start. Thus, these periods get the event times `-k` to `-1`, and
#' the event time 0 is the true start. The cohort is also the period of the
#' true start. A cohort that starts after the last period of the panel is not
#' in a group.
#'
#' To keep the anticipation cells out of the estimate, also remove
#' them from the `target` of the fit. Use [trop_rollout()] to make `W` and
#' `target` from the start dates.
#'
#' The function finds the start of each run from `W`. If a run starts in the
#' first period of the panel, its anticipation periods can be before the panel.
#' Then the event times of the run can be wrong, and the function gives a
#' warning message. If the fit has event times from the argument `event_time`
#' of [trop()], the function uses them. Then `anticipation` must be 0. Refer to
#' the section Event times of [trop()].
#'
#' @section Weights:
#' By default, each cell in a group has the same weight. Thus, a unit with more
#' cells in a group has more weight in the group. With `weights = "unit"`, each
#' unit has the same total weight in each group. The function divides the
#' weight of a unit by its number of cells in the group.
#'
#' The weights are relative. In each group, the function divides the weights
#' by the largest weight of the units of that group. Thus, the estimate of a
#' group does not depend on the weights of other units. In a group, the
#' largest positive weight must be less than 1e200 times the smallest positive
#' weight.
#'
#' With `unit_weights`, each unit has its own weight, for example its number of
#' workers before the treatment. Use weights that the treatment cannot change.
#' With `weights = "cell"`, each cell gets the weight of its unit. With
#' `weights = "unit"`, the weight of a unit is its total weight in each group.
#' A unit with a weight of 0 is not in a group.
#'
#' The bootstrap uses the same weights. A unit that a draw selects two times
#' counts as two units in the draw. The ATT of the fit always gives the same
#' weight to each cell.
#'
#' @section Inference:
#' The bootstrap uses the unit draws of the fit. For each draw, the function
#' calculates the local estimator again with the tuning parameters of the fit.
#' Then it averages the effects of each group. The tuning parameters do not
#' change. The standard error of a group is the standard deviation of its draws.
#'
#' A group that has no treated cell in a draw has no value for that draw. The
#' column `n_draws` gives the number of draws with a value for each group.
#'
#' The draws run in parallel with the `threads` of the control of the fit.
#'
#' By default, the function uses the draws that the fit stored in `fit$boot$rows`.
#' Then the average of all treated cells has the same draws as `fit$boot$draws`.
#' With `n_boot`, the function makes new draws. With `boot_rows`, you give the
#' draws. A group of one fit and a group of another fit have a joint
#' distribution only if they use the same draws.
#'
#' Each call calculates the estimator again for each draw. To get many sets of
#' groups, give them in one list to `by`. To calculate the draws only one time,
#' make the fit with `se = "none"`. Then include `"overall"` in the list, and
#' give the draws with `boot_rows`.
#'
#' The columns `band.low` and `band.high` give a simultaneous band with the
#' sup-t method. The band has a critical value for all groups. The critical
#' value is a quantile of the largest absolute value of the standardized draws
#' over the groups. The band contains all true group effects with the
#' probability `conf_level`, if the bootstrap is correct.
#'
#' The band uses only the draws with a value for all groups. The centers and
#' the standard errors of the band also come from these draws. Thus, if some
#' draws have no value for a group, the band and the interval of that group use
#' different draws. If fewer than 2 draws have a value for all groups, the
#' band is not available.
#'
#' The intervals `conf.low` and `conf.high` are for one group each. The band
#' needs many draws. With only one group, the band columns are `NA`.
#'
#' @section Diagnostics:
#' The attribute `diagnostics` is a data frame with one row for each draw. Its
#' columns count these cells and fits of the draw:
#' * `n_used`: the treated cells with an effect.
#' * `n_skipped`: the cells that the untreated cells do not identify.
#' * `n_nonconverged`: the fits that did not converge.
#' * `n_failed`: the fits that failed. A fit fails if the singular value
#'   decomposition fails or if the result is not finite.
#' * `n_dropped`: the cells that have an effect in the fit, but that the
#'   untreated cells of the draw do not identify. The count does not include
#'   the fits that failed.
#'
#' The function gives a warning message when a draw has such fits or cells.
#' It also gives a warning message for a group with fewer than 2 draws with a
#' value. The print method shows these counts.
#'
#' @section Long runs:
#' With `progress = TRUE`, the function shows the number of draws that are
#' done and an estimate of the time that is left. With `checkpoint`, the
#' function saves the draws to a file after each block of draws. If the
#' calculation stops, call the function again with the same arguments and the
#' same file. The function then calculates only the draws that are not in the
#' file. The results are equal to the results of one run without a stop.
#'
#' The file has a fingerprint of the calculation. The fingerprint uses the
#' panel, the tuning parameters, the solver settings, the draw rows, the groups,
#' the weights, and the version of fasttrop. If the fingerprint of the file is
#' different, the function stops with an error and does not change the file.
#' The number of threads is not in the fingerprint, because it does not change
#' the results.
#'
#' The file is an RDS file with the draws, but not the panel. Keep it in a
#' folder on your computer. The function makes sure that it can write to the
#' folder before it starts the bootstrap.
#'
#' @param fit A `trop` object from [trop()] or [trop_matrix()] with
#'   `method = "local"`.
#' @param by The groups. Refer to the section Groups.
#' @param window An optional vector with 2 event times. The function uses cells
#'   with an event time from `window[1]` to `window[2]`.
#' @param n_boot The number of new bootstrap draws. With `NULL`, the function
#'   uses the draws of the fit, if the fit has draws. With 0, the function does
#'   not calculate a bootstrap.
#' @param boot_rows A matrix of draws from [trop_boot_rows()]. This argument
#'   replaces `n_boot`.
#' @param strata,resample The settings for new draws. Refer to
#'   [trop_boot_rows()]. With `resample = NULL`, the function uses the
#'   setting of the fit.
#' @param conf_level The confidence level. The default is the level of the fit.
#' @param ci Use `"normal"` or `"percentile"`. Refer to [trop()].
#' @param seed An optional integer seed for new draws. The function does not
#'   change the state of the random number generator in the R session.
#' @param anticipation The number of treated periods at the start of each run
#'   that come before the true start. Refer to the section Anticipation. The
#'   default 0 uses the first treated period as the start.
#' @param weights Use `"cell"` or `"unit"`. Refer to the section Weights.
#' @param unit_weights An optional numeric vector with one weight of 0 or more
#'   for each unit. A named vector must have a value for each unit name of the
#'   fit. A vector without names must be in the order of the rows of the panel.
#' @param checkpoint An optional path of a file for the bootstrap draws. Refer
#'   to the section Long runs.
#' @param progress Logical. If `TRUE`, the function shows the progress of the
#'   bootstrap.
#'
#' @return A data frame of class `trop_effects` with one row for each group.
#'   The columns are the group (the name of the column is the value of `by`),
#'   `estimate`, `n_cells`, and `n_units`. With a bootstrap, the columns
#'   `std.error`, `conf.low`, `conf.high`, `band.low`, `band.high`, and
#'   `n_draws` are also in the data frame. The data frame has these attributes:
#'   `draws` (a matrix with one column for each group), `rows` (the rows of the
#'   draws), `critical` (the critical value of the band), `n_complete` (the
#'   number of draws with a value for all groups), `diagnostics`, and `units`
#'   (the names of the units of the fit). When `by` is a list, the result is a
#'   named list of such data frames.
#'
#'   `coef()` gives the estimates with the group names. `vcov()` gives the
#'   covariance matrix of the draws, from the draws with a value for all groups.
#' @seealso [trop_cells()], [trop_contrast()], [trop_boot_rows()],
#'   [trop_stepdown()]
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' fit <- trop(sim, "y", "d", "unit", "time",
#'             lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'             n_boot = 50, seed = 1, control = trop_control(threads = 1))
#' trop_effects(fit, by = "event_time")
#'
#' # Two sets of groups with one bootstrap sweep.
#' both <- trop_effects(fit, by = list(event = "event_time", time = "time"))
#' trop_stepdown(both)
#' @export
trop_effects <- function(fit, by = c("event_time", "cohort", "time", "unit", "overall"),
                         window = NULL, n_boot = NULL, boot_rows = NULL, strata = NULL,
                         resample = NULL, conf_level = fit$conf.level,
                         ci = c("normal", "percentile"), seed = NULL, anticipation = 0L,
                         weights = c("cell", "unit"), unit_weights = NULL,
                         checkpoint = NULL, progress = FALSE) {
  check_local_fit(fit, "trop_effects")
  ci <- match.arg(ci)
  weights <- match.arg(weights)
  check_flag(progress, "progress")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be more than 0 and less than 1.", call. = FALSE)
  }
  P <- fit$panel
  N <- nrow(P$W)
  if (!is.null(strata) && !is.null(boot_rows)) {
    stop("Use `strata` or `boot_rows`, not both. `boot_rows` has the draws, and `strata` makes them.",
         call. = FALSE)
  }
  anticipation <- check_anticipation(anticipation)
  check_event_source(fit, anticipation)
  wt <- effect_weights(weights, unit_weights, fit)
  if (!is.null(checkpoint)) checkpoint <- check_checkpoint_path(checkpoint)
  multi <- is.list(by) && !is.data.frame(by)
  specs <- if (multi) effect_specs(by) else list(by)
  layers <- if (multi) names(specs) else layer_name(by)
  uses_event <- !is.null(window) ||
    any(vapply(specs, function(b) is.character(b) && b[1L] %in% c("event_time", "cohort"), logical(1)))
  if (uses_event) anticipation_warning(P$W, anticipation, P$event)
  groups <- lapply(seq_along(specs), function(k) {
    effect_groups(specs[[k]], W = P$W, tau = fit$tau, window = window, units = fit$units,
                  times = fit$times, anticipation = anticipation, keep = wt$keep,
                  event_given = P$event, layer = layers[k])
  })
  check_weight_range(wt$u, lapply(groups, `[[`, "id"))

  rows <- NULL
  if (!is.null(boot_rows)) {
    rows <- check_boot_rows(boot_rows, N)
  } else if (is.null(n_boot)) {
    rows <- fit$boot$rows
  } else if (n_boot > 0) {
    check_count(n_boot, "n_boot")
    if (!is.null(strata)) strata <- check_strata(strata, N)
    resample <- if (is.null(resample)) {
      if (identical(fit$boot$resample, "controls")) "controls" else "stratified"
    } else {
      match.arg(resample, c("stratified", "controls"))
    }
    treated <- which(rowSums(P$W == 1L & is.finite(P$Y)) > 0L)
    rng <- seed_scope(seed)
    on.exit(rng$restore(), add = TRUE)
    rows <- boot_rows_matrix(treated, setdiff(seq_len(N), treated), as.integer(n_boot), resample,
                             strata)
  }
  if (is.null(rows) && !is.null(checkpoint)) {
    stop("`checkpoint` is only for a bootstrap. Give `n_boot` or `boot_rows`, or use a fit with draws.",
         call. = FALSE)
  }

  sizes <- vapply(groups, function(g) length(g$labels), integer(1))
  boot <- if (is.null(rows)) NULL else {
    effects_boot(fit, lapply(groups, `[[`, "id"), sizes, rows, wt, checkpoint, progress)
  }
  ends <- cumsum(sizes)
  meta <- list(anticipation = anticipation, window = window, weights = wt$label,
               event_given = !is.null(P$event), units = fit$units)
  out <- lapply(seq_along(groups), function(k) {
    d <- if (is.null(boot)) NULL else boot$draws[, (ends[k] - sizes[k] + 1L):ends[k], drop = FALSE]
    effects_table(groups[[k]], fit$tau, d, rows, conf_level, ci, wt, boot$diagnostics, meta)
  })
  if (!is.null(boot)) few_draws_warning(out, if (multi) names(specs) else NULL)
  if (!multi) return(out[[1L]])
  names(out) <- names(specs)
  out
}

layer_name <- function(by) {
  if (is.character(by)) by[1L] else "by"
}

check_event_source <- function(fit, anticipation) {
  if (!is.null(fit$panel$event) && anticipation != 0L) {
    stop("The fit has event times from `event_time`. Then `anticipation` must be 0, because the event times already count the anticipation periods.",
         call. = FALSE)
  }
  invisible(NULL)
}

# A run that starts in the first period can have anticipation periods before
# the panel. Then the event times of the run are not certain.
anticipation_warning <- function(W, anticipation, event_given) {
  if (anticipation == 0L || !is.null(event_given)) return(invisible(NULL))
  n <- sum(W[, 1L] == 1L)
  if (n > 0L) {
    warning(sprintf("%d unit(s) are treated in the first period. Their anticipation periods can be before the panel. Thus, their event times can be wrong. To give the event times, use `event_time` in trop().",
                    n), call. = FALSE)
  }
  invisible(NULL)
}

# A warning for the groups with fewer than 2 draws with a value.
few_draws_warning <- function(out, layers) {
  bad <- unlist(lapply(seq_along(out), function(k) {
    x <- out[[k]]
    lab <- effect_labels(x)[x$n_draws < 2L]
    if (length(lab) == 0L) return(character())
    if (is.null(layers)) lab else paste0(layers[k], ": ", lab)
  }))
  if (length(bad) > 0L) {
    warning(sprintf("%d group(s) have fewer than 2 draws with a value. These groups have no standard error and no interval: %s.",
                    length(bad), paste(utils::head(bad, 6L), collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

check_local_fit <- function(fit, what) {
  if (!inherits(fit, "trop")) stop("`fit` must be a `trop` object.", call. = FALSE)
  if (is.null(fit$panel)) {
    stop("`fit` does not contain the panel. Calculate the fit again with this version of fasttrop.",
         call. = FALSE)
  }
  if (fit$method != "local") {
    stop(sprintf("`%s()` is only for fits with `method = \"local\"`.", what), call. = FALSE)
  }
  invisible(fit)
}

check_anticipation <- function(anticipation) {
  if (!is.numeric(anticipation) || length(anticipation) != 1L || !is.finite(anticipation) ||
      anticipation < 0 || anticipation != round(anticipation) || anticipation > 1e6) {
    stop("`anticipation` must be one whole number from 0 to 1e6.", call. = FALSE)
  }
  as.integer(anticipation)
}

# The weights of the units, from the rule and the optional unit weights. The
# result has the weight of each unit (NULL for 1), the rule, the units with a
# positive weight, and a label for the print method.
effect_weights <- function(weights, unit_weights, fit) {
  N <- nrow(fit$panel$W)
  per_unit <- identical(weights, "unit")
  u <- NULL
  if (!is.null(unit_weights)) {
    if (!is.numeric(unit_weights) || anyNA(unit_weights) || !all(is.finite(unit_weights)) ||
        any(unit_weights < 0)) {
      stop("`unit_weights` must contain finite numbers of 0 or more, and no `NA` values.",
           call. = FALSE)
    }
    nm <- names(unit_weights)
    if (!is.null(nm)) {
      if (is.null(fit$units)) {
        stop("The fit has no unit names. Give `unit_weights` without names, in the order of the rows of the panel.",
             call. = FALSE)
      }
      if (anyDuplicated(nm)) stop("The names of `unit_weights` must be different.", call. = FALSE)
      key <- unit_keys(fit$units, N)
      if (is.null(key)) {
        stop("The unit names of the fit are missing or not different. Give `unit_weights` without names, in the order of the rows of the panel.",
             call. = FALSE)
      }
      miss <- setdiff(key, nm)
      if (length(miss) > 0L) {
        stop(sprintf("`unit_weights` has no value for %d unit(s), for example %s.", length(miss),
                     paste(utils::head(miss, 3L), collapse = ", ")), call. = FALSE)
      }
      unit_weights <- unit_weights[key]
    } else if (length(unit_weights) != N) {
      stop(sprintf("`unit_weights` must have one value for each of the %d units.", N), call. = FALSE)
    }
    if (!any(unit_weights > 0)) stop("All values of `unit_weights` are 0.", call. = FALSE)
    u <- unname(as.double(unit_weights))
  }
  label <- if (is.null(u)) {
    if (per_unit) "each unit has the same total weight in each group" else "each cell has the same weight"
  } else {
    if (per_unit) "the total weight of each unit in each group is its value in `unit_weights`"
    else "the weight of each cell is the value in `unit_weights` of its unit"
  }
  list(u = u, per_unit = per_unit, keep = if (is.null(u)) rep(TRUE, N) else u > 0, label = label)
}

# The weights are relative. In each group, the C code and group_estimates()
# divide the weights by the largest weight of the units of that group. Thus,
# the weights of a group do not depend on the units of other groups. The
# positive weights of a group must have a range that the sums can hold. The
# argument ids is a list of matrices with the group of each cell, or NA.
check_weight_range <- function(w, ids) {
  if (is.null(w)) return(invisible(NULL))
  for (id in ids) {
    on <- !is.na(id)
    if (!any(on)) next
    wi <- w[row(id)[on]]
    gi <- id[on]
    pos <- wi > 0
    if (!any(pos)) next
    lo <- tapply(wi[pos], gi[pos], min)
    hi <- tapply(wi[pos], gi[pos], max)
    if (any(lo / hi < 1e-200)) {
      stop("In a group, the positive values of `unit_weights` are too different. The largest value must be less than 1e200 times the smallest positive value.",
           call. = FALSE)
    }
  }
  invisible(NULL)
}

# The named elements of a list `by`.
effect_specs <- function(by) {
  if (length(by) == 0L) stop("When `by` is a list, it must have at least 1 element.", call. = FALSE)
  nm <- names(by)
  if (is.null(nm)) nm <- rep("", length(by))
  for (k in seq_along(by)) {
    b <- by[[k]]
    ok <- (is.character(b) && length(b) == 1L) || is.matrix(b)
    if (!ok) {
      stop("Each element of the list `by` must be one group name or a matrix.", call. = FALSE)
    }
    if (!nzchar(nm[k])) {
      if (is.matrix(b)) stop("Give a name to each matrix in the list `by`.", call. = FALSE)
      nm[k] <- b
    }
  }
  if (anyDuplicated(nm)) stop("The names of the list `by` must be different.", call. = FALSE)
  names(by) <- nm
  by
}

# Weighted means of the effects in each group.
group_estimates <- function(id, tau, G, wt) {
  use <- is.finite(tau) & !is.na(id)
  gi <- id[use]
  ui <- row(tau)[use]
  w <- if (is.null(wt$u)) rep(1, length(gi)) else wt$u[ui]
  # The weights of a group are divided by its largest weight, as in the C code.
  if (!is.null(wt$u) && length(gi) > 0L) w <- w / stats::ave(w, gi, FUN = max)
  if (wt$per_unit && length(gi) > 0L) w <- w / stats::ave(gi, ui, gi, FUN = length)
  est <- rep(NA_real_, G)
  if (length(gi) > 0L) {
    num <- rowsum(w * tau[use], gi)
    den <- rowsum(w, gi)
    ix <- as.integer(rownames(num))
    est[ix] <- num[, 1L] / den[, 1L]
    bad <- ix[!is.finite(est[ix])]
    if (length(bad) > 0L) {
      # A sum of very large effects can be infinite. Then use scaled effects.
      m <- max(abs(tau[use]))
      num2 <- rowsum(w * (tau[use] / m), gi)
      est[bad] <- m * (num2[match(bad, ix), 1L] / den[match(bad, ix), 1L])
    }
  }
  pairs <- !duplicated(cbind(ui, gi))
  list(estimate = est, n_cells = tabulate(gi, G), n_units = tabulate(gi[pairs], G))
}

# The result for one set of groups.
effects_table <- function(grp, tau, draws, rows, conf_level, ci, wt, diagnostics, meta) {
  G <- length(grp$labels)
  ge <- group_estimates(grp$id, tau, G, wt)
  out <- data.frame(group = grp$labels, estimate = ge$estimate, n_cells = ge$n_cells,
                    n_units = ge$n_units)
  names(out)[1L] <- grp$name
  if (!is.null(draws)) {
    colnames(draws) <- as.character(grp$labels)
    out <- add_boot_columns(out, ge$estimate, draws, rows, conf_level, ci)
    attr(out, "diagnostics") <- diagnostics
  }
  attr(out, "by") <- grp$name
  attr(out, "conf.level") <- conf_level
  attr(out, "anticipation") <- meta$anticipation
  attr(out, "window") <- meta$window
  attr(out, "weights") <- meta$weights
  attr(out, "event_given") <- meta$event_given
  attr(out, "units") <- meta$units
  class(out) <- c("trop_effects", "data.frame")
  out
}

# Standard errors, intervals, the band, and the number of draws for each row.
add_boot_columns <- function(out, est, draws, rows, conf_level, ci) {
  se <- apply(draws, 2L, function(v) {
    v <- v[is.finite(v)]
    if (length(v) >= 2L) stats::sd(v) else NA_real_
  })
  a <- (1 - conf_level) / 2
  out$std.error <- se
  if (ci == "normal") {
    out$conf.low <- est - stats::qnorm(1 - a) * se
    out$conf.high <- est + stats::qnorm(1 - a) * se
  } else {
    q <- apply(draws, 2L, function(v) {
      v <- v[is.finite(v)]
      if (length(v) >= 2L) stats::quantile(v, c(a, 1 - a), type = 7, names = FALSE) else c(NA_real_, NA_real_)
    })
    out$conf.low <- q[1L, ]
    out$conf.high <- q[2L, ]
  }
  band <- if (ncol(draws) > 1L) supt_band(draws, conf_level) else list(crit = NA_real_, se = NA_real_)
  out$band.low <- est - band$crit * band$se
  out$band.high <- est + band$crit * band$se
  out$n_draws <- as.integer(colSums(is.finite(draws)))
  attr(out, "draws") <- draws
  attr(out, "rows") <- rows
  attr(out, "critical") <- band$crit
  attr(out, "n_complete") <- sum(rowSums(!is.finite(draws)) == 0L)
  attr(out, "ci_type") <- ci
  out
}

#' @rdname trop_effects
#' @param x A `trop_effects` object.
#' @param digits The number of significant digits to show.
#' @param row.names,optional Not used.
#' @param ... Not used.
#' @export
print.trop_effects <- function(x, digits = 4L, ...) {
  by <- attr(x, "by")
  cat(sprintf("TROP effects %s\n", group_title(by)))
  for (line in effects_notes(x)) cat("  ", line, "\n", sep = "")
  df <- as.data.frame(x)
  print(df, digits = digits, row.names = FALSE)
  rows <- attr(x, "rows")
  if (!is.null(rows)) {
    B <- nrow(rows)
    nc <- attr(x, "n_complete")
    if (nrow(x) > 1L && is.finite(attr(x, "critical"))) {
      cat(sprintf("Bootstrap: %d draws. The band is simultaneous (sup-t) at the %g%% level. It uses the %d draws with a value for all rows.\n",
                  B, 100 * attr(x, "conf.level"), nc))
      if (nc < B) {
        cat("The intervals use all the draws with a value for each row. Thus, the band and the intervals use different draws.\n")
      }
    } else if (nrow(x) > 1L) {
      cat(sprintf("Bootstrap: %d draws. The band is not available: %d draw(s) have a value for all rows, and the band needs at least 2 draws with a value that changes.\n",
                  B, nc))
    } else {
      cat(sprintf("Bootstrap: %d draws. The band needs 2 or more rows.\n", B))
    }
    lines <- diag_lines(attr(x, "diagnostics"))
    if (length(lines) > 0L) {
      cat(paste0("Bootstrap fits: ", lines, "\n"), sep = "")
      cat("Refer to attr(x, \"diagnostics\") for the counts of each draw.\n")
    }
  }
  invisible(x)
}

# The title of a print: "by event time", or "of the overall estimate".
group_title <- function(by) {
  if (identical(by, "overall")) "of the overall estimate" else paste("by", gsub("_", " ", by, fixed = TRUE))
}

effects_notes <- function(x) {
  out <- character()
  k <- attr(x, "anticipation")
  by <- attr(x, "by")
  if (!is.null(k) && k > 0L && by %in% c("event_time", "cohort")) {
    out <- c(out, sprintf("Event time 0 is the first period after %d anticipation period(s).", k))
  }
  w <- attr(x, "window")
  if (!is.null(w)) out <- c(out, sprintf("Window: event times %s to %s.", format(w[1L]), format(w[2L])))
  if (isTRUE(attr(x, "event_given")) && (by %in% c("event_time", "cohort") || !is.null(w))) {
    out <- c(out, "The event times come from `event_time` of the fit.")
  }
  if (!is.null(attr(x, "weights"))) out <- c(out, sprintf("Weights: %s.", attr(x, "weights")))
  if (!is.null(attr(x, "contrast_note"))) out <- c(out, attr(x, "contrast_note"))
  out
}

#' @rdname trop_effects
#' @export
as.data.frame.trop_effects <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x
  class(out) <- "data.frame"
  for (a in c("draws", "rows", "critical", "by", "conf.level", "n_complete", "diagnostics",
              "anticipation", "window", "weights", "ci_type", "contrast_note", "seconds",
              "event_given", "units", "coefficients")) {
    attr(out, a) <- NULL
  }
  rownames(out) <- NULL
  out
}

#' @rdname trop_effects
#' @param object A `trop_effects` object.
#' @export
coef.trop_effects <- function(object, ...) {
  stats::setNames(object$estimate, effect_labels(object))
}

#' @rdname trop_effects
#' @export
vcov.trop_effects <- function(object, ...) {
  D <- effect_draws(object)
  if (is.null(D)) {
    stop("The object has no bootstrap draws. Calculate it with a bootstrap.", call. = FALSE)
  }
  keep <- rowSums(!is.finite(D)) == 0L
  if (sum(keep) < 2L) {
    stop("Fewer than 2 draws have a value for all rows. The covariance is not available.",
         call. = FALSE)
  }
  V <- stats::cov(D[keep, , drop = FALSE])
  lab <- effect_labels(object)
  dimnames(V) <- list(lab, lab)
  attr(V, "n_draws") <- sum(keep)
  V
}

effect_labels <- function(x) as.character(x[[1L]])

# The draws of a trop_effects object, after a check that the rows of the
# object match the columns of the draws. A subset of the rows keeps the
# attributes, so the draws of other rows would stay in the object.
effect_draws <- function(x) {
  D <- attr(x, "draws")
  if (is.null(D)) return(NULL)
  if (ncol(D) != nrow(x) || !identical(colnames(D), effect_labels(x))) {
    stop("The rows of the `trop_effects` object do not match its draws. Do not remove rows of the object. Use `window` in trop_effects(), or use trop_contrast().",
         call. = FALSE)
  }
  D
}

# Find the group of each treated cell. The result has the group numbers (NA
# for a cell that is in no group) and the labels of the groups. Only cells with
# an effect in `tau` make groups. With `anticipation = k`, the first k periods
# of each run get negative event times, and the cohort is the period after them.
# Units with keep = FALSE are in no group. With event_given, the event times
# come from this matrix, and the cohort is the period with the event time 0.
effect_groups <- function(by, W, tau, window, units, times, anticipation = 0L,
                          keep = rep(TRUE, nrow(W)), event_given = NULL, layer = NULL,
                          allow_empty = FALSE) {
  N <- nrow(W)
  TT <- ncol(W)
  if (is.null(event_given)) {
    start <- run_start(W)
    event <- matrix(NA_integer_, N, TT)
    has <- !is.na(start)
    event[has] <- (col(W) - start)[has] - anticipation
    cohort <- start + anticipation
  } else {
    event <- event_given
    event[W != 1L] <- NA_integer_
    cohort <- col(W) - event
  }
  cohort[!is.na(cohort) & (cohort > TT | cohort < 1L)] <- NA
  if (is.matrix(by)) {
    if (!all(dim(by) == dim(W))) {
      stop("When `by` is a matrix, it must have the dimensions of the panel.", call. = FALSE)
    }
    key <- by
    name <- "group"
  } else {
    name <- by <- match.arg(by, c("event_time", "cohort", "time", "unit", "overall"))
    key <- switch(by, event_time = event, cohort = cohort, time = col(W) * 1L,
                  unit = row(W) * 1L, overall = matrix(1L, N, TT))
  }
  key[W != 1L | !is.finite(tau)] <- NA
  if (!is.null(window)) {
    if (!is.numeric(window) || length(window) != 2L || anyNA(window) || window[1L] > window[2L]) {
      stop("`window` must be 2 numbers, from the first to the last event time.", call. = FALSE)
    }
    key[is.na(event) | event < window[1L] | event > window[2L]] <- NA
  }
  n_before <- sum(!is.na(key))
  key[!keep[row(W)]] <- NA
  lev <- sort(unique(key[!is.na(key)]))
  if (length(lev) == 0L && !allow_empty) {
    where <- if (is.null(layer)) "" else sprintf(" of `%s`", layer)
    if (n_before > 0L) {
      stop(sprintf("All treated cells with an effect in the groups%s have a weight of 0 in `unit_weights`.",
                   where), call. = FALSE)
    }
    stop(sprintf("No treated cell with an effect is in a group%s. Use trop_cells() to count the cells.",
                 where), call. = FALSE)
  }
  labels <- switch(name,
    cohort = , time = if (is.null(times)) lev else times[lev],
    unit = if (is.null(units)) lev else units[lev],
    overall = if (length(lev) > 0L) "overall" else character(),
    lev)
  list(id = matrix(match(key, lev), N, TT), labels = labels, name = name, event = event)
}

# For each treated cell, the first period of the current run of treated
# periods. Other cells have NA.
run_start <- function(W) {
  TT <- ncol(W)
  out <- matrix(NA_integer_, nrow(W), TT)
  for (i in seq_len(nrow(W))) {
    w <- W[i, ]
    st <- which(w == 1L & c(0L, w[-TT]) == 0L)
    if (length(st) > 0L) {
      out[i, w == 1L] <- st[findInterval(which(w == 1L), st)]
    }
  }
  out
}

# Draw the group means. The C code calculates the local estimator again for
# each draw, in parallel, and averages the effects of all sets of groups. With
# a checkpoint or progress messages, the draws run in blocks.
effects_boot <- function(fit, ids, sizes, rows, wt, checkpoint = NULL, progress = FALSE) {
  P <- fit$panel
  N <- nrow(P$W)
  TT <- ncol(P$W)
  offset <- cumsum(c(0L, sizes))[seq_along(ids)]
  gid <- unlist(lapply(seq_along(ids), function(k) {
    g <- as.vector(ids[[k]]) - 1L + offset[k]
    g[is.na(g)] <- -1L
    g
  }), use.names = FALSE)
  storage.mode(gid) <- "integer"
  G <- as.integer(sum(sizes))
  B <- nrow(rows)
  obs <- matrix(as.integer(is.finite(P$Y)), N, TT)
  r0 <- rows - 1L
  storage.mode(r0) <- "integer"
  uw <- if (is.null(wt$u)) NULL else as.double(wt$u)
  ctrl <- unclass(fit$control)
  eligible <- matrix(as.integer(is.finite(fit$tau)), N, TT)
  run <- function(idx) {
    .Call(C_bootstrap, P$Y, P$W, obs, P$X, isTRUE(P$resid), 0L, 0L, fit$lambda,
          r0[idx, , drop = FALSE], ctrl, eligible, gid, G, uw, isTRUE(wt$per_unit))
  }
  t0 <- Sys.time()
  if (is.null(checkpoint) && !progress) {
    bs <- run(seq_len(B))
    state <- list(att = bs$draws, group_draws = bs$group_draws, diag = bs$diagnostics)
  } else {
    key <- if (is.null(checkpoint)) NA_character_ else {
      effects_fingerprint(fit, rows, gid, G, wt)
    }
    state <- if (is.null(checkpoint)) new_boot_state(key, B, G) else checkpoint_open(checkpoint, key, B, G)
    todo <- which(!state$done)
    threads <- as.integer(fit$control$threads)
    block <- max(8L, 4L * threads)
    if (progress) {
      message(sprintf("Bootstrap: %d draws with %d thread(s). %d draw(s) are done.", B, threads,
                      B - length(todo)))
    }
    done_here <- 0L
    for (idx in split(todo, ceiling(seq_along(todo) / block))) {
      bs <- withCallingHandlers(run(idx), error = function(e) {
        if (!is.null(checkpoint)) {
          message(sprintf("The checkpoint file has %d of %d draws. Call the function again with the same arguments to continue.",
                          sum(state$done), B))
        }
      })
      state$att[idx] <- bs$draws
      state$group_draws[idx, ] <- bs$group_draws
      state$diag[idx, ] <- bs$diagnostics
      state$done[idx] <- TRUE
      done_here <- done_here + length(idx)
      if (!is.null(checkpoint)) checkpoint_save(state, checkpoint)
      if (progress) {
        el <- as.double(difftime(Sys.time(), t0, units = "secs"))
        left <- el / done_here * (length(todo) - done_here)
        message(sprintf("Bootstrap draws: %d of %d done. Time: %s. Time left: about %s.",
                        sum(state$done), B, format_seconds(el), format_seconds(left)))
      }
    }
  }
  failed <- sum(!is.finite(state$att))
  if (failed > 0L) {
    warning(sprintf("In %d of %d bootstrap draws, no treated cell has an effect. These draws have no value for any group.",
                    failed, B), call. = FALSE)
  }
  diagnostics <- boot_diagnostics(state$diag)
  boot_warnings(diagnostics)
  list(draws = state$group_draws, diagnostics = diagnostics,
       seconds = as.double(difftime(Sys.time(), t0, units = "secs")))
}

# The sup-t band from the draws with a value for all groups. The centers, the
# standard errors, and the critical value come from these draws. The critical
# value is a quantile of the largest absolute standardized draw over the
# groups. A group with no change in these draws has a standard error of 0.
supt_band <- function(draws, level) {
  none <- list(crit = NA_real_, se = rep(NA_real_, ncol(draws)))
  complete <- rowSums(!is.finite(draws)) == 0L
  if (sum(complete) < 2L) return(none)
  D <- draws[complete, , drop = FALSE]
  se <- apply(D, 2L, stats::sd)
  use <- is.finite(se) & se > 0
  if (!any(use)) return(none)
  Z <- abs(sweep(sweep(D[, use, drop = FALSE], 2L, colMeans(D[, use, drop = FALSE])), 2L,
                 se[use], "/"))
  m <- apply(Z, 1L, max)
  list(crit = unname(stats::quantile(m, level, type = 7)), se = unname(se))
}

#' Romano-Wolf stepdown p-values
#'
#' `trop_stepdown()` adjusts the p-values of many tests for multiple testing.
#' The tests can be for many outcomes, or for many groups of one outcome. The
#' function uses the stepdown method of Romano and Wolf (2005). It controls the
#' probability of one or more false rejections. The method uses the joint
#' bootstrap distribution of all estimates. Thus, all estimates must have the
#' same bootstrap draws.
#'
#' For each estimate, the function calculates the statistic \eqn{t = \hat\theta /
#' se}, where the standard error is the standard deviation of the draws. The
#' null hypothesis is that the effect is 0. The statistics of the draws are the
#' absolute values of the centered draws divided by the standard errors. The
#' function sorts the estimates by \eqn{|t|}. In each step, it compares \eqn{|t|}
#' with the largest statistic of the draws over the estimates that remain.
#'
#' The p-value is (1 + the number of draws with a value that is equal or
#' larger) / (the number of draws + 1). An adjusted p-value is never smaller
#' than the adjusted p-value of an estimate with a larger \eqn{|t|}.
#'
#' The function uses only the draws that have a value for all estimates. It
#' gives a warning message if some draws do not have a value for all
#' estimates. The print method shows the number of draws that it uses. The
#' result is not exact if the number of draws is small.
#'
#' @section Joint draws:
#' Draw b of each object must use the same units. Give the same `boot_rows`
#' from [trop_boot_rows()] to each [trop()] and [trop_effects()] call. Use the
#' same order of the units in each panel.
#'
#' If each object has valid unit names, the function compares the units of
#' each draw. A valid name is not `NA`, and each unit has a different name. A
#' fit of [trop()] always has valid names. A fit of [trop_matrix()] has the row
#' names of `Y`. If no object has unit names, the function compares only the
#' row numbers.
#'
#' If only some objects have valid names, draw b of all objects must have the
#' same row numbers. The objects with valid names must also have the same
#' units. Then the function gives a warning.
#'
#' @param x A `trop` object, a `trop_effects` object with a bootstrap, or a list
#'   of such objects. Use a named list to name the tests. All objects must have
#'   the same bootstrap draws. Refer to the section Joint draws. The result of
#'   [trop_contrast()] is also a `trop_effects` object.
#'
#' @return A data frame of class `trop_stepdown` with one row for each estimate
#'   and with the columns `name`, `estimate`, `std.error`, `t`, `p_raw`, and
#'   `p_adj`. The column `p_raw` is the bootstrap p-value of one test, and
#'   `p_adj` is the adjusted p-value. An estimate with no standard error has no
#'   p-value. The attribute `n_draws` is the number of draws that the function
#'   uses, and `n_boot` is the number of draws.
#' @references Romano, J. P., and Wolf, M. (2005). Stepwise Multiple Testing as
#'   Formalized Data Snooping. *Econometrica*, 73(4), 1237-1282.
#' @seealso [trop_boot_rows()], [trop_effects()], [trop_contrast()]
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' W <- tapply(sim$d, list(sim$unit, sim$time), sum)
#' rows <- trop_boot_rows(W, n_boot = 40, seed = 1)
#' fit <- trop(sim, "y", "d", "unit", "time",
#'             lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'             boot_rows = rows, control = trop_control(threads = 1))
#' trop_stepdown(trop_effects(fit, by = "event_time"))
#' @export
trop_stepdown <- function(x) {
  st <- stack_draws(x, "trop_stepdown")
  res <- romano_wolf(st$estimate, st$draws)
  out <- data.frame(name = st$labels, res, row.names = NULL)
  B <- nrow(st$draws)
  n_used <- attr(res, "n_draws")
  if (n_used < B) {
    warning(sprintf("%d of %d draws do not have a value for all estimates. The stepdown uses the other %d draws.",
                    B - n_used, B, n_used), call. = FALSE)
  }
  attr(out, "n_draws") <- n_used
  attr(out, "n_boot") <- B
  class(out) <- c("trop_stepdown", "data.frame")
  out
}

#' @rdname trop_stepdown
#' @param digits The number of significant digits to show.
#' @param ... Not used.
#' @export
print.trop_stepdown <- function(x, digits = 4L, ...) {
  cat("Romano-Wolf stepdown p-values\n")
  df <- x
  class(df) <- "data.frame"
  attr(df, "n_draws") <- NULL
  attr(df, "n_boot") <- NULL
  print(df, digits = digits, row.names = FALSE)
  cat(sprintf("Draws: %d of %d have a value for all estimates. The p-values use these draws.\n",
              attr(x, "n_draws"), attr(x, "n_boot")))
  invisible(x)
}

# Put the estimates and the draws of fits and trop_effects objects together.
# Draw b of all objects must use the same units. The labels name the estimates.
stack_draws <- function(x, what) {
  if (inherits(x, "trop") || inherits(x, "trop_effects")) x <- list(x)
  if (!is.list(x) || length(x) == 0L) {
    stop("`x` must be a `trop` object, a `trop_effects` object, or a list of these objects.",
         call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm)) nm <- rep("", length(x))
  est <- numeric(0)
  lab <- character(0)
  draws <- NULL
  rows <- NULL
  units <- NULL
  level <- NULL
  keyed <- TRUE
  weak <- character()
  named <- FALSE
  all_rows <- all_ids <- vector("list", length(x))
  for (k in seq_along(x)) {
    e <- x[[k]]
    own <- if (nzchar(nm[k])) nm[k] else paste0("test", k)
    if (inherits(e, "trop")) {
      r <- e$boot$rows
      d <- if (is.null(e$boot$draws)) NULL else matrix(e$boot$draws, ncol = 1L)
      u <- e$units
      est_k <- e$estimate
      lab_k <- own
      level <- level %||% e$conf.level
    } else if (inherits(e, "trop_effects")) {
      r <- attr(e, "rows")
      d <- effect_draws(e)
      u <- attr(e, "units")
      est_k <- e$estimate
      g <- effect_labels(e)
      lab_k <- if (!nzchar(nm[k]) && length(x) == 1L) g else paste0(own, ": ", g)
      level <- level %||% attr(e, "conf.level")
    } else {
      stop("Each element of `x` must be a `trop` object or a `trop_effects` object.", call. = FALSE)
    }
    if (is.null(r) || is.null(d)) {
      stop("Each object must have a bootstrap. Calculate it with `se = \"bootstrap\"` or with `n_boot`.",
           call. = FALSE)
    }
    keys <- unit_keys(u, ncol(r))
    if (!is.null(u)) named <- TRUE
    if (is.null(keys)) {
      keyed <- FALSE
      weak <- c(weak, own)
    }
    ids_k <- if (is.null(keys)) NULL else matrix(keys[r], nrow(r), ncol(r))
    all_rows[[k]] <- r
    all_ids[k] <- list(ids_k)
    if (k == 1L) {
      rows <- r
      units <- u
    }
    est <- c(est, est_k)
    lab <- c(lab, lab_k)
    draws <- cbind(draws, d)
  }
  # Without names for all objects, the function cannot compare the units. If
  # no object has names, this is the documented use, and there is no warning.
  if (length(x) > 1L && length(weak) > 0L && named) {
    warning(sprintf("These objects have no valid unit names: %s. A valid name is not NA, and each unit has a different name. Thus, draw b of each object must have the same row numbers. Make sure that each panel has the same order of the units.",
                    paste(weak, collapse = ", ")), call. = FALSE)
  }
  # The objects with valid names must have the same units in each draw. If an
  # object has no valid names, all objects must also have the same row numbers.
  same <- function(a, b) identical(dim(a), dim(b)) && all(a == b)
  has_ids <- !vapply(all_ids, is.null, logical(1))
  ok <- TRUE
  if (sum(has_ids) > 1L) {
    ref <- all_ids[[which(has_ids)[1L]]]
    ok <- all(vapply(all_ids[has_ids], same, logical(1), b = ref))
  }
  if (ok && !all(has_ids)) {
    ok <- all(vapply(all_rows, same, logical(1), b = all_rows[[1L]]))
  }
  if (!ok) {
    stop("All objects must use the same bootstrap draws. Draw b of each object must have the same units. Give the same `boot_rows` from `trop_boot_rows()` to each call, and use the same order of the units in each panel.",
         call. = FALSE)
  }
  if (anyDuplicated(lab)) {
    stop("The estimates must have different names. Give a named list.", call. = FALSE)
  }
  colnames(draws) <- lab
  # A result of these objects has the unit names only if each object has valid
  # names. Otherwise, a later comparison can use names that do not apply.
  list(estimate = est, draws = draws, labels = lab, rows = rows,
       units = if (keyed) units else NULL, conf_level = level)
}

# The unit names as text for the comparison of draws. The names are valid only
# if each unit has a name that is not NA and that is different from the names
# of the other units. A number becomes text that gives the same number again.
unit_keys <- function(u, N) {
  if (is.null(u) || length(u) != N) return(NULL)
  # The text of a number must not depend on the decimal mark of the session.
  op <- options(OutDec = ".")
  on.exit(options(op), add = TRUE)
  if (is.factor(u)) u <- as.character(u)
  if (is.double(u) && !is.object(u)) {
    if (anyNA(u) || any(!is.finite(u))) return(NULL)
    u[u == 0] <- 0
    k <- as.character(u)
    back <- suppressWarnings(as.numeric(k))
    long <- is.na(back) | back != u
    k[long] <- sprintf("%.17g", u[long])
  } else {
    # Integers, text, and classes such as integer64 have exact text.
    k <- as.character(u)
  }
  if (anyNA(k) || anyDuplicated(k)) return(NULL)
  k
}

`%||%` <- function(a, b) if (is.null(a)) b else a

romano_wolf <- function(est, draws) {
  K <- length(est)
  D <- draws[rowSums(!is.finite(draws)) == 0L, , drop = FALSE]
  B <- nrow(D)
  if (B < 2L) stop("The objects have fewer than 2 bootstrap draws with a value for all estimates.", call. = FALSE)
  se <- apply(D, 2L, stats::sd)
  Z <- abs(sweep(sweep(D, 2L, colMeans(D)), 2L, se, "/"))
  tstat <- est / se
  p_raw <- p_adj <- rep(NA_real_, K)
  usable <- which(is.finite(se) & se > 0 & is.finite(est))
  if (length(usable) > 0L) {
    ord <- usable[order(abs(tstat[usable]), decreasing = TRUE)]
    prev <- 0
    for (j in seq_along(ord)) {
      rest <- ord[j:length(ord)]
      mx <- if (length(rest) == 1L) Z[, rest] else apply(Z[, rest, drop = FALSE], 1L, max)
      prev <- max(prev, (1 + sum(mx >= abs(tstat[ord[j]]))) / (B + 1))
      p_adj[ord[j]] <- prev
    }
    p_raw[usable] <- vapply(usable, function(k) (1 + sum(Z[, k] >= abs(tstat[k]))) / (B + 1), numeric(1))
  }
  out <- data.frame(estimate = est, std.error = unname(se), t = unname(tstat), p_raw = p_raw,
                    p_adj = p_adj)
  attr(out, "n_draws") <- B
  out
}
