# Tests for trop_rollout(), trop_cells(), weights, checkpoints, contrasts, and
# the inputs for a sensitivity analysis.

l_w <- c(time = 0.3, unit = 0.5, nn = 1)

# A staggered panel in long format with start months and 1 anticipation month.
rollout_panel <- function(seed = 21, N = 18, TT = 10) {
  set.seed(seed)
  units <- sprintf("u%02d", seq_len(N))
  starts <- c(rep(c(5, 7, 9), each = 3), rep(NA, N - 9))
  panel <- expand.grid(unit = units, month = seq_len(TT), stringsAsFactors = FALSE)
  panel$start <- starts[match(panel$unit, units)]
  e <- panel$month - panel$start
  panel$y <- rnorm(N)[match(panel$unit, units)] + cumsum(rnorm(TT, 0, 0.3))[panel$month] +
    rnorm(nrow(panel), 0, 0.2) + ifelse(!is.na(e) & e >= 0, 1 + 0.3 * e, 0)
  panel
}

test_that("trop_rollout makes the treatment, the target, and the event time", {
  panel <- data.frame(unit = rep(c("b", "a", "c"), each = 5), month = rep(c(2001, 2003, 2002, 2005, 2004), 3))
  panel$start <- c(b = 2003, a = 2004, c = NA)[panel$unit]
  r <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 1, window = c(0, 1))
  expect_named(r, c("treated", "target", "event_time"))
  b <- r[panel$unit == "b", ][order(panel$month[panel$unit == "b"]), ]
  expect_equal(b$event_time, -2:2)
  expect_equal(b$treated, c(0L, 1L, 1L, 1L, 1L))
  expect_equal(b$target, c(FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_true(all(is.na(r$event_time[panel$unit == "c"])))
  expect_true(all(r$treated[panel$unit == "c"] == 0L & !r$target[panel$unit == "c"]))
  lead <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 2, window = c(-2, Inf))
  expect_equal(lead$target, lead$treated == 1L)
  expect_error(trop_rollout(panel$unit, panel$month, replace(panel$start, 1, 2002)), "one value for each unit")
  expect_error(trop_rollout(panel$unit, panel$month, ifelse(is.na(panel$start), NA, 1999)), "value of `time`")
  expect_error(trop_rollout(c(panel$unit, "b"), c(panel$month, 2001), c(panel$start, 2003)), "two or more rows")
  expect_error(trop_rollout(panel$unit, panel$month[-1], panel$start), "same length")
  expect_error(trop_rollout(panel$unit, panel$month, panel$start, anticipation = -1), "anticipation")
  expect_error(trop_rollout(panel$unit, panel$month, panel$start, window = c(2, 1)), "window")
})

test_that("trop_cells counts every treated cell once and labels anticipation", {
  panel <- rollout_panel()
  r <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 1, window = c(0, 2))
  panel <- cbind(panel, r)
  panel$y[panel$unit == "u02" & panel$month == 6] <- NA
  fit <- suppressWarnings(trop(panel, "y", "treated", "unit", "month", target = "target",
                               lambda = l_w, se = "none", control = small_control()))
  cells <- trop_cells(fit, anticipation = 1)
  expect_s3_class(cells, "trop_cells")
  expect_equal(cells$event_time, -1:5)
  expect_equal(cells$n_excluded[cells$event_time == -1], 9L)
  expect_equal(cells$n_cells[cells$event_time == -1], 0L)
  expect_equal(sum(cells$n_cells), fit$n_treated_cells - fit$n_skipped)
  expect_equal(sum(cells$share), 1)
  tot <- attr(cells, "totals")
  expect_equal(tot[["with_outcome"]], sum(tot[c("in_estimate", "excluded", "skipped", "zero_weight",
                                                "outside_window", "no_group")]))
  expect_equal(tot[["no_outcome"]], 1)
  expect_equal(tot[["excluded"]], fit$n_excluded_cells)
  expect_output(print(cells), "Event time 0 is the first period after 1 anticipation period")
  expect_output(print(cells), "1 with no outcome")

  win <- trop_cells(fit, anticipation = 1, window = c(0, 1))
  expect_equal(win$event_time, 0:1)
  expect_gt(attr(win, "totals")[["outside_window"]], 0)
  ev <- trop_effects(fit, by = "event_time", anticipation = 1, n_boot = 0)
  expect_equal(ev$n_cells, cells$n_cells[cells$n_cells > 0])
  co <- trop_cells(fit, by = "cohort", anticipation = 1)
  expect_equal(co$cohort, c(5, 7, 9))
  expect_error(trop_cells(fit, by = list("cohort")), "not a list")
})

test_that("weights give the means of a hand calculation in the estimate and in the draws", {
  panel <- rollout_panel(seed = 23)
  panel$treated <- as.integer(!is.na(panel$start) & panel$month >= panel$start)
  fit <- suppressWarnings(trop(panel, "y", "treated", "unit", "month", lambda = l_w,
                               se = "none", control = small_control()))
  u <- setNames(c(3, 1, 2, 5, 0, 1, 2, 2, 4, rep(1, 9)), sprintf("u%02d", 1:18))
  tau <- fit$tau
  hand <- function(tau, w, per_unit) {
    i <- row(tau)[is.finite(tau)]
    v <- tau[is.finite(tau)]
    cw <- w[i]
    if (per_unit) cw <- cw / ave(i, i, FUN = length)
    sum(cw * v) / sum(cw)
  }
  for (per_unit in c(FALSE, TRUE)) {
    rule <- if (per_unit) "unit" else "cell"
    ov <- trop_effects(fit, by = "overall", weights = rule, unit_weights = u, n_boot = 0)
    expect_equal(ov$estimate, hand(tau, unname(u), per_unit), tolerance = 1e-12)
    expect_equal(ov$n_units, 8L)
    eq <- trop_effects(fit, by = "overall", weights = rule, n_boot = 0)
    expect_equal(eq$estimate, hand(tau, rep(1, 18), per_unit), tolerance = 1e-12)
  }
  expect_false(isTRUE(all.equal(
    trop_effects(fit, by = "overall", weights = "unit", n_boot = 0)$estimate, fit$estimate)))

  # The draws: each resampled row is a unit with the weight of its original unit.
  rows <- trop_boot_rows(fit$panel$W, 4, seed = 5)
  rows[1, 1:3] <- c(1L, 1L, 4L)
  ev <- trop_effects(fit, by = "overall", weights = "unit", unit_weights = u, boot_rows = rows)
  for (b in seq_len(nrow(rows))) {
    rb <- rows[b, ]
    inp <- fasttrop:::as_panel_inputs(fit$panel$Y[rb, ], fit$panel$W[rb, ])
    tb <- .Call(fasttrop:::C_local_estimate, inp$Y, inp$W, inp$obs, NULL, fit$lambda,
                unclass(fit$control))$tau
    expect_equal(unname(attr(ev, "draws")[b, 1L]), hand(tb, unname(u)[rb], TRUE), tolerance = 1e-12)
  }
  expect_output(print(ev), "Weights: the total weight of each unit")

  expect_error(trop_effects(fit, unit_weights = u[-1]), "no value for 1 unit")
  expect_error(trop_effects(fit, unit_weights = unname(u)[-1]), "one value for each of the 18 units")
  expect_error(trop_effects(fit, unit_weights = replace(u, 1, -1)), "0 or more")
  expect_error(trop_effects(fit, unit_weights = u * 0), "All values")
  expect_error(trop_effects(fit, weights = "worker"), "should be one of")
})

test_that("a checkpoint resumes a bootstrap with the same results", {
  panel <- rollout_panel(seed = 25)
  panel$treated <- as.integer(!is.na(panel$start) & panel$month >= panel$start)
  fit <- suppressWarnings(trop(panel, "y", "treated", "unit", "month", lambda = l_w,
                               se = "none", control = small_control()))
  rows <- trop_boot_rows(fit$panel$W, 20, seed = 9)
  by <- list(event = "event_time", cohort = "cohort")
  plain <- trop_effects(fit, by = by, boot_rows = rows)
  ck <- tempfile(fileext = ".rds")
  on.exit(unlink(ck), add = TRUE)
  msgs <- capture_messages(trop_effects(fit, by = by, boot_rows = rows, checkpoint = ck,
                                        progress = TRUE))
  expect_match(msgs, "Bootstrap: 20 draws with 2 thread", all = FALSE)
  expect_match(msgs, "Bootstrap draws: 20 of 20 done", all = FALSE)
  st <- readRDS(ck)
  expect_true(all(st$done))
  # A file of a run that stopped: the draws that are not done have no values.
  st$done[c(2, 11:20)] <- FALSE
  st$att[c(2, 11:20)] <- NA
  st$group_draws[c(2, 11:20), ] <- NA
  st$diag[c(2, 11:20), ] <- 0L
  st$check <- payload_hash(st)
  saveRDS(st, ck)
  msgs <- capture_messages(res <- trop_effects(fit, by = by, boot_rows = rows, checkpoint = ck,
                                               progress = TRUE))
  expect_match(msgs, "9 draw\\(s\\) are done", all = FALSE)
  expect_match(msgs, "Bootstrap draws: 17 of 20 done", all = FALSE)
  for (k in names(by)) {
    expect_identical(attr(res[[k]], "draws"), attr(plain[[k]], "draws"))
    expect_identical(as.data.frame(res[[k]]), as.data.frame(plain[[k]]))
  }
  # A file of another calculation, with other groups, gives an error and does not change.
  before <- readRDS(ck)
  expect_error(trop_effects(fit, by = by, window = c(0, 1), boot_rows = rows, checkpoint = ck),
               "different calculation")
  expect_identical(readRDS(ck), before)
  bad <- tempfile(fileext = ".rds")
  writeLines("text", bad)
  expect_error(trop_effects(fit, boot_rows = rows, checkpoint = bad), "cannot read the checkpoint")
  saveRDS(list(a = 1), bad)
  expect_error(trop_effects(fit, boot_rows = rows, checkpoint = bad), "not a fasttrop checkpoint")
  unlink(bad)
  expect_error(trop_effects(fit, boot_rows = rows, checkpoint = file.path(tempdir(), "no-dir", "x.rds")),
               "folder of the checkpoint file does not exist")
  expect_error(trop_effects(fit, checkpoint = ck), "only for a bootstrap")
  # The number of threads is not in the fingerprint.
  fit1 <- fit
  fit1$control <- trop_control(threads = 1L)
  expect_identical(attr(trop_effects(fit1, by = by, boot_rows = rows, checkpoint = ck)$event, "draws"),
                   attr(plain$event, "draws"))
})

test_that("trop_contrast, coef, and vcov use the joint draws", {
  p <- small_panel(seed = 27, n_units = 20, n_treated = 6)
  fit <- trop_matrix(p$Y, p$W, lambda = l_w, n_boot = 25, seed = 2, control = small_control())
  ev <- trop_effects(fit, by = "event_time")
  D <- attr(ev, "draws")
  expect_equal(coef(ev), setNames(ev$estimate, c("0", "1", "2")))
  V <- vcov(ev)
  expect_equal(unname(V), unname(stats::cov(D)), ignore_attr = TRUE)
  expect_equal(attr(V, "n_draws"), 25L)
  ct <- trop_contrast(ev, list(mean = c("0" = 1, "1" = 1, "2" = 1) / 3, slope = c("2" = 1, "0" = -1)))
  expect_s3_class(ct, "trop_effects")
  expect_equal(ct$estimate, c(mean(ev$estimate), ev$estimate[3] - ev$estimate[1]), tolerance = 1e-12)
  expect_equal(unname(attr(ct, "draws")[, "slope"]), unname(D[, 3] - D[, 1]), tolerance = 1e-12)
  expect_equal(ct$std.error[2], stats::sd(D[, 3] - D[, 1]), tolerance = 1e-12)
  expect_equal(nrow(trop_stepdown(list(ev = ev, c = ct))), 5L)
  both <- trop_effects(fit, by = list(ov = "overall", ev = "event_time"))
  d2 <- trop_contrast(both, c("ev: 0" = 1, "ov: overall" = -1))
  expect_equal(d2$estimate, ev$estimate[1] - fit$estimate, tolerance = 1e-12)
  expect_error(trop_contrast(ev, c("7" = 1)), "not estimates: 7")
  expect_error(trop_contrast(ev, c("0" = 0)), "only coefficients of 0")
  expect_error(trop_contrast(ev, list(c("0" = 1))), "must have a name")
  expect_error(trop_contrast(ev, c(1, 2)), "names and finite values")
  # A subset of the rows does not match the draws.
  expect_error(vcov(ev[1:2, ]), "do not match its draws")
  expect_error(trop_stepdown(ev[-1, ]), "do not match its draws")
})

test_that("trop_honest_inputs gives estimates relative to the reference", {
  sim <- trop_simulate(n_units = 30, n_periods = 12, n_treated = 8, n_post = 3, seed = 4)
  start <- ave(ifelse(sim$d == 1, sim$time, NA), sim$unit,
               FUN = function(v) if (all(is.na(v))) NA else min(v, na.rm = TRUE))
  r <- trop_rollout(sim$unit, sim$time, start, anticipation = 3, window = c(-3, Inf))
  sim$placebo <- r$treated
  sim$tg <- r$target
  fit <- trop(sim, "y", "placebo", "unit", "time", target = "tg", lambda = l_w, n_boot = 30,
              seed = 1, control = small_control())
  ev <- trop_effects(fit, by = "event_time", anticipation = 3)
  expect_equal(ev$event_time, -3:2)
  inp <- trop_honest_inputs(ev, reference = -1)
  D <- attr(ev, "draws")
  expect_equal(inp$numPrePeriods, 2L)
  expect_equal(inp$numPostPeriods, 3L)
  expect_equal(inp$event_time, c(-3L, -2L, 0L, 1L, 2L))
  expect_equal(unname(inp$betahat), ev$estimate[-3] - ev$estimate[3], tolerance = 1e-12)
  expect_equal(unname(inp$sigma), unname(stats::cov(D[, -3] - D[, 3])), tolerance = 1e-12)
  ref2 <- trop_honest_inputs(ev, reference = -2)
  expect_equal(ref2$numPrePeriods, 1L)
  expect_equal(ref2$numPostPeriods, 4L)
  expect_error(trop_honest_inputs(trop_effects(fit, by = "cohort")), "by = \"event_time\"")
  expect_error(trop_honest_inputs(ev, reference = 0), "less than 0")
  expect_error(trop_honest_inputs(ev, reference = -5), "not in `x`")
  expect_error(trop_honest_inputs(ev, reference = -3), "before the reference")
  expect_error(trop_honest_inputs(trop_effects(fit, by = "event_time", anticipation = 3, n_boot = 0)),
               "bootstrap draws")
  holes <- fit
  holes$tau[, 8] <- NA
  expect_error(trop_honest_inputs(trop_effects(holes, by = "event_time", anticipation = 3)),
               "no gaps")
})
