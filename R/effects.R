#' Effects by event time, cohort, period, or unit
#'
#' `trop_effects()` averages the effects of the treated cells of a [trop()] fit
#' in groups. A group is, for example, all treated cells with the same number
#' of periods since the start of the treatment. The function can also give
#' bootstrap standard errors, intervals, and simultaneous bands for all groups.
#'
#' The function uses the effects `tau` of the fit. Thus, the target cells of the
#' fit are the only cells in the groups. A fit that has an `NA` effect for a
#' cell does not give a value for that cell.
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
#'
#' The argument `window` limits the cells to cells with an event time in the
#' window. For example, `window = c(0, 3)` uses the first 4 treated periods.
#'
#' @section Inference:
#' The bootstrap uses the unit draws of the fit. For each draw, the function
#' calculates the local estimator again with the tuning parameters of the fit.
#' Then it averages the effects of each group. The tuning parameters do not
#' change. The standard error of a group is the standard deviation of its draws.
#' A group that has no treated cell in a draw has no value for that draw.
#'
#' By default, the function uses the draws that the fit stored in `fit$boot$rows`.
#' Then the average of all treated cells has the same draws as `fit$boot$draws`.
#' With `n_boot`, the function makes new draws. With `boot_rows`, you give the
#' draws. A group of one fit and a group of another fit have a joint
#' distribution only if they use the same draws.
#'
#' The columns `band.low` and `band.high` give a simultaneous band with the
#' sup-t method. The band has a critical value for all groups. The critical
#' value is a quantile of the largest absolute value of the standardized draws
#' over the groups. The band contains all true group effects with the
#' probability `conf_level`, if the bootstrap is correct.
#'
#' The intervals `conf.low` and `conf.high` are for one group each. The band
#' needs many draws. With only one group, the band columns are `NA`.
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
#'
#' @return A data frame of class `trop_effects` with one row for each group.
#'   The columns are the group (the name of the column is the value of `by`),
#'   `estimate`, and `n_cells`. With a bootstrap, the columns `std.error`,
#'   `conf.low`, `conf.high`, `band.low`, and `band.high` are also in the data
#'   frame. The data frame has these attributes: `draws` (a matrix with one
#'   column for each group), `rows` (the rows of the draws), and `critical` (the
#'   critical value of the band).
#' @seealso [trop_boot_rows()], [trop_stepdown()]
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' fit <- trop(sim, "y", "d", "unit", "time",
#'             lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'             n_boot = 50, seed = 1, control = trop_control(threads = 1))
#' trop_effects(fit, by = "event_time")
#' @export
trop_effects <- function(fit, by = c("event_time", "cohort", "time", "unit", "overall"),
                         window = NULL, n_boot = NULL, boot_rows = NULL, strata = NULL,
                         resample = NULL, conf_level = fit$conf.level,
                         ci = c("normal", "percentile"), seed = NULL) {
  if (!inherits(fit, "trop")) stop("`fit` must be a `trop` object.", call. = FALSE)
  if (is.null(fit$panel)) {
    stop("`fit` does not contain the panel. Calculate the fit again with this version of fasttrop.",
         call. = FALSE)
  }
  if (fit$method != "local") {
    stop("`trop_effects()` is only for fits with `method = \"local\"`.", call. = FALSE)
  }
  ci <- match.arg(ci)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be more than 0 and less than 1.", call. = FALSE)
  }
  P <- fit$panel
  N <- nrow(P$W)
  if (!is.null(strata) && !is.null(boot_rows)) {
    stop("Use `strata` or `boot_rows`, not both. `boot_rows` has the draws, and `strata` makes them.",
         call. = FALSE)
  }
  grp <- effect_groups(by, P$W, window, fit$units, fit$times)
  G <- length(grp$labels)

  tau <- fit$tau
  use <- is.finite(tau) & !is.na(grp$id)
  est <- rep(NA_real_, G)
  n_cells <- tabulate(grp$id[use], G)
  s <- rowsum(tau[use], grp$id[use])
  ix <- as.integer(rownames(s))
  est[ix] <- s[, 1L] / n_cells[ix]

  out <- data.frame(group = grp$labels, estimate = est, n_cells = n_cells)
  names(out)[1L] <- grp$name

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

  attr(out, "by") <- grp$name
  attr(out, "conf.level") <- conf_level
  if (!is.null(rows)) {
    draws <- effects_boot(fit, grp$id, G, rows)
    colnames(draws) <- as.character(grp$labels)
    se <- apply(draws, 2L, function(v) {
      v <- v[is.finite(v)]
      if (length(v) >= 2L) stats::sd(v) else NA_real_
    })
    a <- (1 - conf_level) / 2
    if (ci == "normal") {
      out$std.error <- se
      out$conf.low <- est - stats::qnorm(1 - a) * se
      out$conf.high <- est + stats::qnorm(1 - a) * se
    } else {
      q <- apply(draws, 2L, function(v) {
        v <- v[is.finite(v)]
        if (length(v) >= 2L) stats::quantile(v, c(a, 1 - a), type = 7, names = FALSE) else c(NA_real_, NA_real_)
      })
      out$std.error <- se
      out$conf.low <- q[1L, ]
      out$conf.high <- q[2L, ]
    }
    crit <- if (G > 1L) supt_critical(draws, conf_level) else NA_real_
    out$band.low <- est - crit * se
    out$band.high <- est + crit * se
    attr(out, "draws") <- draws
    attr(out, "rows") <- rows
    attr(out, "critical") <- crit
  }
  class(out) <- c("trop_effects", "data.frame")
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
  cat(sprintf("TROP effects by %s\n", gsub("_", " ", by, fixed = TRUE)))
  df <- as.data.frame(x)
  print(df, digits = digits, row.names = FALSE)
  if (!is.null(attr(x, "rows"))) {
    cat(sprintf("Bootstrap: %d draws. The band is simultaneous (sup-t) at the %g%% level.\n",
                nrow(attr(x, "rows")), 100 * attr(x, "conf.level")))
  }
  invisible(x)
}

#' @rdname trop_effects
#' @export
as.data.frame.trop_effects <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x
  class(out) <- "data.frame"
  for (a in c("draws", "rows", "critical", "by", "conf.level")) attr(out, a) <- NULL
  rownames(out) <- NULL
  out
}

# Find the group of each treated cell. The result has the group numbers (NA
# for a cell that is in no group) and the labels of the groups.
effect_groups <- function(by, W, window, units, times) {
  N <- nrow(W)
  TT <- ncol(W)
  start <- run_start(W)
  event <- matrix(NA_integer_, N, TT)
  has <- !is.na(start)
  event[has] <- (col(W) - start)[has]
  if (is.matrix(by)) {
    if (!all(dim(by) == dim(W))) {
      stop("When `by` is a matrix, it must have the dimensions of the panel.", call. = FALSE)
    }
    key <- by
    name <- "group"
  } else {
    name <- by <- match.arg(by, c("event_time", "cohort", "time", "unit", "overall"))
    key <- switch(by, event_time = event, cohort = start, time = col(W) * 1L, unit = row(W) * 1L,
                  overall = matrix(1L, N, TT))
  }
  key[W != 1L] <- NA
  if (!is.null(window)) {
    if (!is.numeric(window) || length(window) != 2L || anyNA(window) || window[1L] > window[2L]) {
      stop("`window` must be 2 numbers, from the first to the last event time.", call. = FALSE)
    }
    key[is.na(event) | event < window[1L] | event > window[2L]] <- NA
  }
  lev <- sort(unique(key[!is.na(key)]))
  if (length(lev) == 0L) stop("No treated cell is in a group.", call. = FALSE)
  labels <- switch(name,
    cohort = , time = if (is.null(times)) lev else times[lev],
    unit = if (is.null(units)) lev else units[lev],
    overall = "overall",
    lev)
  list(id = matrix(match(key, lev), N, TT), labels = labels, name = name)
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

# Draw the group means. Each draw calculates the local estimator again.
effects_boot <- function(fit, id, G, rows) {
  P <- fit$panel
  ctrl <- unclass(fit$control)
  B <- nrow(rows)
  draws <- matrix(NA_real_, B, G)
  failed <- 0L
  for (b in seq_len(B)) {
    rb <- rows[b, ]
    Yb <- P$Y[rb, , drop = FALSE]
    Wb <- P$W[rb, , drop = FALSE]
    ob <- matrix(as.integer(is.finite(Yb)), nrow(Yb), ncol(Yb))
    Xb <- if (is.null(P$X)) NULL else P$X[rb, , , drop = FALSE]
    est <- tryCatch({
      if (isTRUE(P$resid) && !is.null(Xb)) {
        phi <- .Call(C_residualize, Yb, Wb, ob, Xb, ctrl)
        Yb <- subtract_xphi(Yb, Xb, phi)
        Xb <- NULL
      }
      .Call(C_local_estimate, Yb, Wb, ob, Xb, fit$lambda, ctrl)
    }, error = function(e) NULL)
    if (is.null(est)) {
      failed <- failed + 1L
      next
    }
    gb <- id[rb, , drop = FALSE]
    ok <- is.finite(est$tau) & !is.na(gb)
    if (!any(ok)) next
    n <- tabulate(gb[ok], G)
    s <- rowsum(est$tau[ok], gb[ok])
    ix <- as.integer(rownames(s))
    draws[b, ix] <- s[, 1L] / n[ix]
  }
  if (failed > 0L) {
    warning(failed, " of ", B, " bootstrap draws failed. The standard errors do not use these draws.",
            call. = FALSE)
  }
  draws
}

# The critical value of the sup-t band: a quantile of the largest absolute
# standardized draw over the groups.
supt_critical <- function(draws, level) {
  se <- apply(draws, 2L, function(v) {
    v <- v[is.finite(v)]
    if (length(v) >= 2L) stats::sd(v) else NA_real_
  })
  use <- is.finite(se) & se > 0
  if (!any(use)) return(NA_real_)
  Z <- abs(sweep(sweep(draws[, use, drop = FALSE], 2L, colMeans(draws[, use, drop = FALSE], na.rm = TRUE)),
                 2L, se[use], "/"))
  keep <- rowSums(is.finite(Z)) > 0L
  if (sum(keep) < 2L) return(NA_real_)
  m <- apply(Z[keep, , drop = FALSE], 1L, max, na.rm = TRUE)
  unname(stats::quantile(m, level, type = 7))
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
#' The p-value is (1 + the number of draws with a larger value) / (the number of
#' draws + 1). An adjusted p-value is never smaller than the adjusted p-value of
#' an estimate with a larger \eqn{|t|}.
#'
#' The function uses only the draws that have a value for all estimates. The
#' result is not exact if the number of draws is small.
#'
#' @param x A `trop` object, a `trop_effects` object with a bootstrap, or a list
#'   of such objects. Use a named list to name the tests. All objects must have
#'   the same bootstrap rows. To make them, give the same `boot_rows` from
#'   [trop_boot_rows()] to each [trop()] and [trop_effects()] call.
#'
#' @return A data frame with one row for each estimate and with the columns
#'   `name`, `estimate`, `std.error`, `t`, `p_raw`, and `p_adj`. The column
#'   `p_raw` is the bootstrap p-value of one test, and `p_adj` is the adjusted
#'   p-value. An estimate with no standard error has no p-value.
#' @references Romano, J. P., and Wolf, M. (2005). Stepwise Multiple Testing as
#'   Formalized Data Snooping. *Econometrica*, 73(4), 1237-1282.
#' @seealso [trop_boot_rows()], [trop_effects()]
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
  for (k in seq_along(x)) {
    e <- x[[k]]
    own <- if (nzchar(nm[k])) nm[k] else paste0("test", k)
    if (inherits(e, "trop")) {
      r <- e$boot$rows
      d <- if (is.null(e$boot$draws)) NULL else matrix(e$boot$draws, ncol = 1L)
      est_k <- e$estimate
      lab_k <- own
    } else if (inherits(e, "trop_effects")) {
      r <- attr(e, "rows")
      d <- attr(e, "draws")
      est_k <- e$estimate
      g <- as.character(e[[1L]])
      lab_k <- if (!nzchar(nm[k]) && length(x) == 1L) g else paste0(own, ": ", g)
    } else {
      stop("Each element of `x` must be a `trop` object or a `trop_effects` object.", call. = FALSE)
    }
    if (is.null(r) || is.null(d)) {
      stop("Each object must have a bootstrap. Calculate it with `se = \"bootstrap\"` or with `n_boot`.",
           call. = FALSE)
    }
    if (!is.null(rows) && !(all(dim(rows) == dim(r)) && all(rows == r))) {
      stop("All objects must use the same bootstrap rows. Give the same `boot_rows` from `trop_boot_rows()` to each call.",
           call. = FALSE)
    }
    rows <- r
    est <- c(est, est_k)
    lab <- c(lab, lab_k)
    draws <- cbind(draws, d)
  }
  res <- romano_wolf(est, draws)
  data.frame(name = lab, res, row.names = NULL)
}

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
  data.frame(estimate = est, std.error = se, t = tstat, p_raw = p_raw, p_adj = p_adj)
}
