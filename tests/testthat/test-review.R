# Regression tests for the findings of the review of version 0.4.0. Each test
# reproduces one finding with a small panel.

l_r <- c(time = 0.3, unit = 0.5, nn = 1)

# A staggered panel in long format with a start for each unit.
review_panel <- function(seed = 31, N = 16, TT = 10) {
  set.seed(seed)
  units <- sprintf("u%02d", seq_len(N))
  starts <- c(rep(c(4, 6, 8), each = 3), rep(NA, N - 9))
  d <- expand.grid(unit = units, time = seq_len(TT), stringsAsFactors = FALSE)
  d$start <- starts[match(d$unit, units)]
  e <- d$time - d$start
  d$y <- rnorm(N)[match(d$unit, units)] + rnorm(TT)[d$time] + rnorm(nrow(d), 0, 0.2) +
    ifelse(!is.na(e) & e >= 0, 2, 0)
  d
}

test_that("a missing row does not change the event times of a run", {
  d <- review_panel()
  r <- trop_rollout(d$unit, d$time, d$start, anticipation = 2, window = c(0, Inf))
  d <- cbind(d, r)
  # Unit u04 starts in period 6. Remove its rows for periods 4 and 7.
  gone <- d$unit == "u04" & d$time %in% c(4, 7)
  dm <- d[!gone, ]
  w <- collect_warnings(fit <- trop(dm, "y", "treated", "unit", "time", target = "target",
                                    lambda = l_r, se = "none", control = small_control()))
  expect_false(any(grepl("treatment stops", w)))
  expect_true(any(grepl("no rows for the periods just before the first treated row", w)))
  expect_equal(fit$design, "staggered")
  # The interior gap (period 7) does not stop the run.
  expect_equal(fit$panel$W[which(fit$units == "u04"), ], c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1))

  # With event_time, the event times and the cohorts are exactly the ones of
  # trop_rollout(), and there is no warning about the start.
  w2 <- collect_warnings(fit2 <- trop(dm, "y", "treated", "unit", "time", target = "target",
                                      event_time = "event_time", lambda = l_r, se = "none",
                                      control = small_control()))
  expect_false(any(grepl("just before the first treated row", w2)))
  ev <- trop_effects(fit2, by = "event_time", n_boot = 0)
  expect_equal(ev$event_time, 0:6)
  cells <- trop_cells(fit2, by = "cohort")
  expect_equal(cells$cohort, c(4, 6, 8))
  expect_equal(cells$n_units, c(3L, 3L, 3L))
  # The row of u04 for event time -2 is missing. Thus, 8 cells have event time
  # -2 and 9 cells have event time -1.
  ce <- trop_cells(fit2, by = "event_time")
  expect_equal(ce$n_excluded[ce$event_time == -2], 8L)
  expect_equal(ce$n_excluded[ce$event_time == -1], 9L)
  expect_output(print(ev), "The event times come from `event_time` of the fit")
  expect_error(trop_effects(fit2, by = "event_time", anticipation = 2), "must be 0")
  expect_error(trop_cells(fit2, anticipation = 1), "must be 0")
})

test_that("event_time is checked, and a run in the first period gives a warning", {
  d <- review_panel()
  r <- trop_rollout(d$unit, d$time, d$start)
  d <- cbind(d, r)
  d$bad <- d$event_time
  d$bad[d$treated == 1L][1L] <- NA
  expect_error(trop(d, "y", "treated", "unit", "time", event_time = "bad", lambda = l_r,
                    se = "none", control = small_control()), "no value for 1 treated cell")
  d$bad <- d$event_time + 0.5
  expect_error(trop(d, "y", "treated", "unit", "time", event_time = "bad", lambda = l_r,
                    se = "none", control = small_control()), "whole numbers")
  expect_error(trop(d, "y", "treated", "unit", "time", event_time = "nope", lambda = l_r,
                    se = "none", control = small_control()), "name of one column")
  W <- tapply(d$treated, list(d$unit, d$time), sum)
  Y <- tapply(d$y, list(d$unit, d$time), sum)
  expect_error(trop_matrix(Y, W, event_time = matrix(0, 2, 2), lambda = l_r, se = "none"),
               "dimensions of `Y`")

  # Unit 1 has two runs. The true start of the first run is period 1, and the
  # true start of the second run is period 7. The other treated units start in
  # period 7. Each run has 1 anticipation period. Thus, the anticipation period
  # of the first run of unit 1 is before the panel.
  set.seed(57)
  N <- 12
  TT <- 8
  Wg <- matrix(0L, N, TT)
  Wg[1, c(1:2, 6:8)] <- 1L
  Wg[2:4, 6:8] <- 1L
  E <- matrix(NA_real_, N, TT)
  E[1, c(1:2, 6:8)] <- c(0, 1, -1, 0, 1)
  E[2:4, 6:8] <- rep(c(-1, 0, 1), each = 3)
  Yg <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + Wg +
    matrix(rnorm(N * TT, 0, 0.1), N, TT)
  twfe <- c(time = 0, unit = 0, nn = Inf)
  fw <- suppressWarnings(trop_matrix(Yg, Wg, lambda = twfe, se = "none", control = small_control()))
  expect_warning(ev_w <- trop_effects(fw, by = "event_time", anticipation = 1, n_boot = 0),
                 "1 unit\\(s\\) are treated in the first period")
  expect_equal(ev_w$event_time, -1:1)
  expect_equal(ev_w$n_cells, c(5L, 5L, 4L))
  fe <- suppressWarnings(trop_matrix(Yg, Wg, lambda = twfe, se = "none", event_time = E,
                                     control = small_control()))
  ev_e <- trop_effects(fe, by = "event_time", n_boot = 0)
  expect_equal(ev_e$n_cells, c(4L, 5L, 5L))
  expect_equal(ev_e$estimate[ev_e$event_time == 0],
               mean(fe$tau[cbind(c(1, 1, 2, 3, 4), c(1, 7, 7, 7, 7))]))
  co <- trop_cells(fe, by = "cohort")
  expect_equal(co$cohort, c(1L, 7L))
  expect_equal(co$n_cells, c(2L, 12L))
})

test_that("joint inference needs the same units in each draw, not the same row numbers", {
  p <- small_panel(seed = 33, n_units = 16, n_treated = 5)
  Y <- p$Y
  W <- p$W
  rownames(Y) <- rownames(W) <- sprintf("u%02d", seq_len(nrow(Y)))
  rows <- trop_boot_rows(W, 12, seed = 4)
  a <- trop_matrix(Y, W, lambda = l_r, boot_rows = rows, control = small_control())
  # The same panel with the units in another order.
  perm <- c(rev(which(rowSums(W) > 0)), rev(which(rowSums(W) == 0)))
  b_bad <- suppressWarnings(trop_matrix(Y[perm, ], W[perm, ], lambda = l_r, boot_rows = rows,
                                        control = small_control()))
  expect_error(trop_contrast(list(a = a, b = b_bad), c(a = 1, b = -1)), "same units")
  expect_error(trop_stepdown(list(a = a, b = b_bad)), "same units")
  # The rows that select the same units in the new order.
  rows_b <- matrix(match(rownames(Y)[rows], rownames(Y)[perm]), nrow(rows))
  b_ok <- trop_matrix(Y[perm, ], W[perm, ], lambda = l_r, boot_rows = rows_b,
                      control = small_control())
  ct <- trop_contrast(list(a = a, b = b_ok), c(a = 1, b = -1))
  expect_equal(ct$estimate, 0, tolerance = 1e-8)
  expect_lt(ct$std.error, 1e-8)
  ea <- trop_effects(a, by = "overall")
  eb <- trop_effects(b_ok, by = "overall")
  expect_equal(nrow(trop_stepdown(list(a = ea, b = eb))), 2L)
})

test_that("a missing estimate does not change another contrast", {
  p <- small_panel(seed = 35, n_units = 16, n_treated = 5)
  rows <- trop_boot_rows(p$W, 10, seed = 3)
  good <- trop_matrix(p$Y, p$W, lambda = l_r, boot_rows = rows, control = small_control())
  bad <- good
  bad$estimate <- NaN
  ct <- trop_contrast(list(good = good, bad = bad), list(good_only = c(good = 1), bad_only = c(bad = 1),
                                                          both = c(good = 1, bad = 1)))
  expect_equal(ct$estimate, c(good$estimate, NA, NA))
  expect_equal(unname(attr(ct, "draws")[, "good_only"]), good$boot$draws)
})

test_that("the global method skips cells that the draw does not identify", {
  set.seed(37)
  N <- 12
  TT <- 8
  W <- matrix(0L, N, TT)
  W[1:4, 7:8] <- 1L
  Y <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + 2 * W
  # Only unit 12 is an untreated unit with an outcome in period 7.
  Y[5:11, 7] <- NA
  twfe <- c(time = 0, unit = 0, nn = Inf)
  rows <- rbind(c(1:4, 5:11, 12), c(1:4, 5:11, 11))
  w <- collect_warnings(gl <- trop_matrix(Y, W, method = "global", lambda = twfe,
                                          boot_rows = rows, control = small_control()))
  expect_equal(gl$estimate, 2, tolerance = 1e-8)
  expect_equal(gl$boot$draws[1], 2, tolerance = 1e-8)
  d <- gl$boot$diagnostics
  expect_equal(d$n_skipped, c(0L, 4L))
  expect_equal(d$n_dropped, c(0L, 4L))
  expect_equal(d$n_used, c(8L, 4L))
  expect_equal(gl$boot$draws[2], 2, tolerance = 1e-8)
  expect_true(any(grepl("do not identify 4 treated cell", w)))
  # The point estimate also skips a cell that the data do not identify.
  Y2 <- Y
  Y2[12, 7] <- NA
  w2 <- collect_warnings(g2 <- trop_matrix(Y2, W, method = "global", lambda = twfe, se = "none",
                                           control = small_control()))
  expect_equal(g2$n_skipped, 4L)
  expect_equal(g2$estimate, 2, tolerance = 1e-8)
  expect_true(any(grepl("did not use 4 treated cell", w2)))
})

test_that("n_dropped counts each cell, and a failure does not cancel a loss", {
  set.seed(39)
  N <- 12
  TT <- 6
  W <- matrix(0L, N, TT)
  W[1:4, 5:6] <- 1L
  Y <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + W
  Y[5:11, 6] <- NA
  Y[1, 4] <- NA
  twfe <- c(time = 0, unit = 0, nn = Inf)
  fit <- trop_matrix(Y, W, lambda = twfe, se = "none", control = small_control())
  rows <- rbind(1:12, c(1:11, 11), c(1, 1, 1, 1, 5:11, 11))
  expect_warning(bs <- trop_effects(fit, by = "overall", boot_rows = rows),
                 "do not identify 8 treated cell")
  d <- attr(bs, "diagnostics")
  expect_equal(d$n_dropped, c(0L, 4L, 4L))
  expect_equal(d$n_skipped, d$n_dropped)
  expect_true(all(d$n_failed == 0L))

  # A failed fit in the same draw does not cancel the loss of a cell. Bit 4 of
  # the fault test makes all SVDs fail in the second draw only.
  ln <- c(time = 0, unit = 0, nn = 1)
  ctl <- small_control()
  ctl$fault <- 7L
  fit2 <- trop_matrix(Y, W, lambda = ln, se = "none", control = ctl)
  expect_equal(sum(is.finite(fit2$tau)), 8L)
  rows2 <- rbind(1:12, c(1:11, 11), c(1:11, 11))
  w <- collect_warnings(b2 <- trop_effects(fit2, by = "overall", boot_rows = rows2))
  d2 <- attr(b2, "diagnostics")
  expect_equal(d2$n_skipped, c(0L, 4L, 4L))
  expect_equal(d2$n_dropped, c(0L, 4L, 4L))
  expect_equal(d2$n_failed, c(0L, 4L, 0L))
  expect_true(any(grepl("fit\\(s\\) failed", w)))
})

test_that("the band uses the complete draws and says when it is not available", {
  s <- local({
    set.seed(41)
    N <- 12
    TT <- 6
    W <- matrix(0L, N, TT)
    W[1:2, 5:6] <- 1L
    W[3:4, 6] <- 1L
    Y <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + W +
      matrix(rnorm(N * TT, 0, 0.1), N, TT)
    list(Y = Y, W = W)
  })
  fit <- suppressWarnings(trop_matrix(s$Y, s$W, lambda = c(time = 0, unit = 0, nn = Inf),
                                      se = "none", control = small_control()))
  # In the first draws, units 1 and 2 are not in the draw. In the other draws,
  # units 3 and 4 are not in the draw. Thus, no draw has a value for all cohorts.
  rows <- rbind(t(replicate(5, c(3, 4, 3, 4, 5:12))), t(replicate(5, c(1, 2, 1, 2, 5:12))))
  w <- collect_warnings(co <- trop_effects(fit, by = "cohort", boot_rows = rows))
  expect_true(any(grepl("draws with a value", w)) || all(co$n_draws >= 2L))
  expect_equal(attr(co, "n_complete"), 0L)
  expect_true(all(is.na(co$band.low)))
  expect_output(print(co), "The band is not available")
})

test_that("a damaged checkpoint file gives an error", {
  p <- small_panel(seed = 43)
  fit <- trop_matrix(p$Y, p$W, lambda = l_r, se = "none", control = small_control())
  rows <- trop_boot_rows(p$W, 10, seed = 2)
  ck <- tempfile(fileext = ".rds")
  on.exit(unlink(ck), add = TRUE)
  trop_effects(fit, by = "overall", boot_rows = rows, checkpoint = ck)
  st <- readRDS(ck)
  bad <- st
  bad$done <- bad$done[1:4]
  saveRDS(bad, ck)
  expect_error(trop_effects(fit, by = "overall", boot_rows = rows, checkpoint = ck), "is damaged")
  bad <- st
  bad$done[3] <- FALSE
  saveRDS(bad, ck)
  expect_error(trop_effects(fit, by = "overall", boot_rows = rows, checkpoint = ck), "is damaged")
  bad <- st
  bad$diag <- bad$diag[, 1:4]
  saveRDS(bad, ck)
  expect_error(trop_effects(fit, by = "overall", boot_rows = rows, checkpoint = ck), "is damaged")
})

test_that("the cross-validation prefers values with a score for more cells", {
  # The outcomes of the units are far apart. With a large unit decay, the unit
  # weights are 0, and many cells have no score. The sum over fewer cells must
  # not win.
  set.seed(45)
  N <- 14
  TT <- 8
  W <- matrix(0L, N, TT)
  W[1:4, 7:8] <- 1L
  Y <- outer(seq_len(N) * 40, rep(1, TT)) + outer(rep(1, N), rnorm(TT)) +
    matrix(rnorm(N * TT), N, TT) + W
  g <- trop_grid(time = 0, unit = c(0, 600), nn = Inf)
  w <- collect_warnings(fit <- trop_matrix(Y, W, grid = g, se = "none", control = small_control()))
  tr <- fit$cv$trace
  big <- tr[tr$unit == 600, ]
  expect_true(all(big$n_valid < max(tr$n_valid)))
  expect_equal(fit$lambda[["unit"]], 0)
  expect_equal(fit$cv$n_valid_chosen, max(tr$n_valid))
  expect_output(print(summary(fit)), sprintf("Tuning score: uses %d of %d cells", max(tr$n_valid),
                                             fit$cv$n_cells))
})

test_that("a grid with no usable score gives an error", {
  p <- small_panel(seed = 47)
  ctl <- small_control()
  ctl$fault <- 3L
  g <- trop_grid(time = c(0, 0.5), unit = 0, nn = c(0.5, 1))
  expect_error(trop_matrix(p$Y, p$W, grid = g, cv_cells = 20, se = "none", seed = 1, control = ctl),
               "no tuning values with a usable score")
  expect_error(trop_matrix(p$Y, p$W, method = "global", grid = g, n_placebo = 10, se = "none",
                           seed = 1, control = ctl), "no tuning values with a usable score")
})

test_that("weights are relative, and the messages name an empty layer", {
  d <- review_panel(seed = 49)
  d$treated <- as.integer(!is.na(d$start) & d$time >= d$start)
  fit <- suppressWarnings(trop(d, "y", "treated", "unit", "time", lambda = l_r, se = "none",
                               control = small_control()))
  units <- as.character(fit$units)
  u <- setNames(rep(1, length(units)), units)
  base <- trop_effects(fit, by = "cohort", weights = "unit", unit_weights = u, n_boot = 0)
  huge <- trop_effects(fit, by = "cohort", weights = "unit", unit_weights = u * 1e308, n_boot = 0)
  expect_equal(huge$estimate, base$estimate)
  rows <- trop_boot_rows(fit$panel$W, 6, seed = 2)
  hb <- trop_effects(fit, by = "overall", unit_weights = u * 1e308, boot_rows = rows)
  expect_true(all(hb$n_draws == 6L))

  z <- u
  z[c("u01", "u02", "u03")] <- 0
  expect_error(trop_effects(fit, by = list(ev = "event_time", co = matrix(ifelse(fit$panel$W == 1 &
                 row(fit$panel$W) <= 3, 1L, NA_integer_), nrow(fit$panel$W))), unit_weights = z, n_boot = 0),
               "weight of 0")
  expect_error(trop_effects(fit, by = list(ev = "event_time", late = "cohort"), window = c(50, 60),
                            n_boot = 0), "in a group of `ev`")
  cells <- trop_cells(fit, by = "event_time", window = c(50, 60))
  expect_equal(nrow(cells), 0L)
  expect_equal(attr(cells, "totals")[["outside_window"]], attr(cells, "totals")[["with_outcome"]])
  expect_output(print(cells), "No treated cell is in a group")
})

test_that("a group with fewer than 2 draws gives a warning", {
  d <- review_panel(seed = 51)
  d$treated <- as.integer(!is.na(d$start) & d$time >= d$start)
  fit <- suppressWarnings(trop(d, "y", "treated", "unit", "time", lambda = l_r, se = "none",
                               control = small_control()))
  units <- as.character(fit$units)
  u <- setNames(rep(1, length(units)), units)
  u[c("u01", "u02", "u03")] <- 0
  u["u04"] <- 1
  # The draws select only units with a weight of 0 among the treated units.
  rows <- t(replicate(4, c(1, 1, 1, 1, 2, 3, 1, 2, 3, 10:16)))
  w <- collect_warnings(e <- trop_effects(fit, by = "overall", unit_weights = u,
                                          boot_rows = rows))
  expect_equal(e$n_draws, 0L)
  expect_true(any(grepl("fewer than 2 draws with a value", w)))
  w2 <- collect_warnings(trop_effects(fit, by = "overall",
                                      boot_rows = t(replicate(3, c(10:16, 10:16, 10, 11)))))
  expect_true(any(grepl("no treated cell has an effect", w2)))
})

test_that("a folder that is not writable gives an error before the bootstrap", {
  p <- small_panel(seed = 53)
  fit <- trop_matrix(p$Y, p$W, lambda = l_r, se = "none", control = small_control())
  rows <- trop_boot_rows(p$W, 4, seed = 1)
  bad <- file.path(tempdir(), "no-such-folder", "x.rds")
  expect_error(trop_effects(fit, by = "overall", boot_rows = rows, checkpoint = bad),
               "folder of the checkpoint file does not exist")
  expect_error(trop_rollout(1:4, 1:4, rep(2, 4), anticipation = Inf), "whole number")
  expect_error(trop_rollout(1:4, 1:4, rep(2, 4), anticipation = 2^31), "whole number")
})

test_that("trop_honest_inputs gives a target vector for event time 0", {
  sim <- trop_simulate(n_units = 30, n_periods = 12, n_treated = 8, n_post = 3, seed = 55)
  start <- ave(ifelse(sim$d == 1, sim$time, NA), sim$unit,
               FUN = function(v) if (all(is.na(v))) NA else min(v, na.rm = TRUE))
  r <- trop_rollout(sim$unit, sim$time, start, anticipation = 3, window = c(-3, Inf))
  sim <- cbind(sim, r)
  fit <- trop(sim, "y", "treated", "unit", "time", target = "target", event_time = "event_time",
              lambda = l_r, n_boot = 20, seed = 1, control = small_control())
  ev <- trop_effects(fit, by = "event_time")
  inp <- trop_honest_inputs(ev, reference = -2)
  expect_equal(names(inp$l_vec), c("-1", "0", "1", "2"))
  expect_equal(unname(inp$l_vec), c(0, 1, 0, 0))
  avg <- trop_honest_inputs(ev, reference = -2, event_times = 0:2)
  expect_equal(unname(avg$l_vec), c(0, 1, 1, 1) / 3)
  expect_error(trop_honest_inputs(ev, reference = -2, event_times = -3), "after the reference")
  expect_error(trop_honest_inputs(ev, reference = -2, event_times = c(0, 0)), "different whole")
})

# Round 2 of the review.

test_that("draws pair by unit names only if the names are valid", {
  expect_null(unit_keys(c("a", "a", "b"), 3))
  expect_null(unit_keys(c("a", NA, "b"), 3))
  expect_null(unit_keys(c("a", "b"), 3))
  expect_equal(unit_keys(c(3, 1, 2), 3), c("3", "1", "2"))
  expect_equal(unit_keys(factor(c("x", "y")), 2), c("x", "y"))
  # Two different numbers have the same short text. Their names are different.
  k <- unit_keys(c(1, 1 + .Machine$double.eps, 3), 3)
  expect_equal(anyDuplicated(k), 0L)
  expect_equal(as.numeric(k), c(1, 1 + .Machine$double.eps, 3))

  p <- small_panel(seed = 33, n_units = 16, n_treated = 5)
  rows <- trop_boot_rows(p$W, 12, seed = 4)
  fit_y <- function(Y, W = p$W) {
    trop_matrix(Y, W, lambda = l_r, boot_rows = rows, control = small_control())
  }
  Y <- p$Y
  rownames(Y) <- sprintf("u%02d", seq_len(nrow(Y)))
  # Two treated rows have the same name. A swap of these rows does not change
  # the names. Thus, the function compares only the row numbers, with a warning.
  tr <- which(rowSums(p$W) > 0)
  Yd <- Y
  rownames(Yd)[tr[1:2]] <- "dup"
  perm <- seq_len(nrow(Y))
  perm[tr[1:2]] <- tr[2:1]
  a <- fit_y(Yd)
  b <- fit_y(Yd[perm, ], p$W[perm, ])
  expect_warning(st <- trop_stepdown(list(a = a, b = b)), "no valid unit names")
  expect_equal(attr(st, "n_draws"), 12L)
  # A missing name gives the same warning, not an error.
  Yna <- Y
  rownames(Yna)[1] <- NA
  na_fit <- fit_y(Yna)
  expect_warning(trop_stepdown(list(a = na_fit, b = na_fit)), "no valid unit names")
  # A combination with an object without names has no names.
  anon <- fit_y(p$Y)
  named <- fit_y(Y)
  expect_warning(ct <- trop_contrast(list(a = anon, b = named), c(a = 1)),
                 "no valid unit names: a")
  expect_null(attr(ct, "units"))
  # Without names for any object, the comparison of row numbers is the
  # documented use, and there is no warning.
  expect_no_warning(trop_stepdown(list(a = anon, b = fit_y(p$Y))))
  ct2 <- trop_contrast(list(a = named, b = named), c(a = 1, b = -1))
  expect_equal(attr(ct2, "units"), rownames(Y))
  expect_equal(ct2$std.error, 0)
})

test_that("event times must increase by 1 in a run of treated periods", {
  set.seed(61)
  N <- 7
  TT <- 8
  W <- matrix(0L, N, TT)
  W[1:3, 4:7] <- 1L
  Y <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + 2 * W
  E <- matrix(NA_real_, N, TT)
  E[1:3, 4:7] <- rep(0:3, each = 3)
  twfe <- c(time = 0, unit = 0, nn = Inf)
  fm <- function(W, E) {
    suppressWarnings(trop_matrix(Y, W, lambda = twfe, se = "none", event_time = E,
                                 control = small_control()))
  }
  expect_s3_class(fm(W, E), "trop")
  E_bad <- E
  E_bad[1, 4:7] <- c(0, -1, 2, 3)
  expect_error(fm(W, E_bad), "do not increase by 1")
  # A new run after an untreated period can start a new count.
  W2 <- W
  W2[1, 6] <- 0L
  E2 <- E
  E2[1, 6] <- NA
  E2[1, 7] <- 0
  co <- trop_cells(fm(W2, E2), by = "cohort")
  expect_equal(co$cohort, c(4L, 7L))
  expect_equal(co$n_cells, c(10L, 1L))
  expect_error(fm(W, replace(E, 1, 1e7)), "from -1e6 to 1e6")
})

test_that("the global cross-validation counts the placebo cells, not the sets", {
  # With a very large unit decay, the weights of the units with other
  # outcomes are 0. Then the untreated cells do not identify some placebo
  # cells, and the RMSE of the sets uses fewer cells.
  N <- 13L
  TT <- 7L
  W <- matrix(0L, N, TT)
  W[1:5, 5:7] <- 1L
  u <- c(rep(0, 9), rep(10, 4))
  Y <- outer(u, 0:(TT - 1) / 4, "+") + 2 * W
  Y[6:13, 5:6] <- Y[6:13, 5:6] + rep(c(-0.3, -0.1, 0.1, 0.3), 2)
  Y[6:13, 7] <- NA_real_
  Y[c(6, 10), 7] <- u[c(6, 10)] + 1.5 + c(100, -100)
  g <- trop_grid(time = 0, unit = c(0, 1000), nn = Inf)
  q <- trop_matrix(Y, W, method = "global", grid = g, n_placebo = 40L, seed = 31,
                   se = "none", control = tight_control())
  tr <- q$cv$trace
  expect_equal(q$cv$n_cells, 90L)
  expect_equal(tr$n_valid[tr$unit == 0], 90L)
  expect_equal(tr$n_valid[tr$unit == 1000], 80L)
  expect_lt(tr$score[tr$unit == 1000], tr$score[tr$unit == 0])
  expect_equal(q$lambda[["unit"]], 0)
  expect_output(print(summary(q)), "uses 90 of 90 placebo cells")
})

test_that("a global draw uses the same treated units as the fit", {
  # Unit 6 has W = 1 but no treated outcome. It is not a treated unit of the
  # fit. Thus, it must not change the weights of a draw.
  Y <- outer(c(rep(0, 5), 60, 0, 0, 1, 1, 1, 1), 0:6 / 2, "+")
  Y[1:5, 5:7] <- Y[1:5, 5:7] + 2
  Y[7:8, 5:7] <- Y[7:8, 5:7] + 4
  Y[6, 5:7] <- NA_real_
  W <- matrix(0L, 12, 7)
  W[1:6, 5:7] <- 1L
  a <- suppressWarnings(trop_matrix(Y, W, method = "global", lambda = c(time = 0, unit = 1, nn = Inf),
                                    boot_rows = rbind(1:12, 1:12), control = tight_control()))
  expect_equal(a$boot$draws, rep(a$estimate, 2))
})

test_that("the checkpoint does not change other files, and it finds changed results", {
  p <- small_panel(seed = 47)
  fit <- trop_matrix(p$Y, p$W, lambda = l_r, se = "none", control = small_control())
  rows <- trop_boot_rows(p$W, 10, seed = 5)
  d <- tempfile("ck-")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ck <- file.path(d, "draws.rds")
  other <- paste0(ck, ".tmp")
  writeLines("other", other)
  run <- function(by = "event_time") trop_effects(fit, by = by, boot_rows = rows, checkpoint = ck)
  ev <- run()
  expect_equal(readLines(other), "other")
  # In this block design, the groups by period have the same cells as the
  # groups by event time. Thus, the draws of the file are correct for them.
  expect_identical(unname(attr(run("time"), "draws")), unname(attr(ev, "draws")))
  expect_error(run("overall"), "different calculation")
  expect_equal(readLines(other), "other")
  expect_setequal(list.files(d, all.files = TRUE, no.. = TRUE), c("draws.rds", "draws.rds.tmp"))

  full <- readRDS(ck)
  # A correct state with 2 draws that are not done. The function calculates
  # these draws again and gets the same results.
  part <- full
  part$done[9:10] <- FALSE
  part$att[9:10] <- NA_real_
  part$group_draws[9:10, ] <- NA_real_
  part$diag[9:10, ] <- 0L
  part$check <- payload_hash(part)
  saveRDS(part, ck)
  expect_identical(attr(run(), "draws"), attr(ev, "draws"))
  # Draws that are not done, with the flag of done draws.
  bad <- part
  bad$done[9:10] <- TRUE
  bad$check <- payload_hash(bad)
  saveRDS(bad, ck)
  expect_error(run(), "is damaged")
  # An ATT that is not finite in a draw that used treated cells.
  bad <- full
  bad$att[1] <- Inf
  bad$check <- payload_hash(bad)
  saveRDS(bad, ck)
  expect_error(run(), "is damaged")
  # A change of a value without a new hash.
  bad <- full
  bad$att[1] <- bad$att[1] + 1
  saveRDS(bad, ck)
  expect_error(run(), "is damaged")
})

test_that("a large weight of a unit outside of the groups does not change the estimates", {
  p <- small_panel(seed = 35, n_units = 16, n_treated = 5)
  rows <- trop_boot_rows(p$W, 10, seed = 3)
  fit <- trop_matrix(p$Y, p$W, lambda = l_r, boot_rows = rows, control = small_control())
  tr <- rowSums(p$W) > 0
  for (rule in c("cell", "unit")) {
    ref <- trop_effects(fit, by = "overall", weights = rule, unit_weights = as.double(tr))
    ext <- trop_effects(fit, by = "overall", weights = rule, unit_weights = ifelse(tr, 1e-300, 1e300))
    expect_identical(ext$estimate, ref$estimate)
    expect_identical(attr(ext, "draws"), attr(ref, "draws"))
    tiny <- ifelse(tr, .Machine$double.xmin * .Machine$double.eps, 1)
    ut <- trop_effects(fit, by = "overall", weights = rule, unit_weights = tiny)
    expect_identical(ut$estimate, ref$estimate)
    expect_equal(trop_cells(fit, by = "overall", weights = rule, unit_weights = tiny)$share, 1)
  }
  w <- as.double(tr)
  w[which(tr)[1]] <- 1e-250
  expect_error(trop_effects(fit, by = "overall", unit_weights = w), "too different")
})

test_that("the summary separates draws without an estimate from failed fits", {
  p <- small_panel(seed = 49, n_units = 14, n_treated = 5)
  tr <- which(rowSums(p$W) > 0)
  ctrl <- which(rowSums(p$W) == 0)
  donor <- seq_len(nrow(p$W))
  donor[tr] <- ctrl[1]
  rows <- rbind(seq_len(nrow(p$W)), donor, seq_len(nrow(p$W)))
  w <- collect_warnings(fit <- trop_matrix(p$Y, p$W, lambda = l_r, boot_rows = rows,
                                           control = small_control()))
  expect_true(any(grepl("no treated cell has an effect", w)))
  expect_equal(fit$boot$diagnostics$n_used[2], 0L)
  expect_equal(fit$boot$diagnostics$n_failed[2], 0L)
  expect_output(print(summary(fit)), "3 draws \\(custom\\), 1 without an estimate")
})

# Round 3 of the review.

test_that("unit names keep the usual text, and only exact text where necessary", {
  k <- unit_keys(c(1e15, 1 + .Machine$double.eps, 2), 3)
  expect_equal(k[c(1, 3)], c("1e+15", "2"))
  expect_equal(as.numeric(k[2]), 1 + .Machine$double.eps)
  expect_equal(unit_keys(c(-0, 1), 2), c("0", "1"))
  expect_equal(unit_keys(c(5L, 7L), 2), c("5", "7"))
  long <- file.path(tempdir(), paste0(strrep("a", 150), ".rds"))
  expect_lt(nchar(basename(own_temp_file(long, ".probe"))), 30)
})

test_that("named weights use the exact unit names of the fit", {
  p <- small_panel(seed = 35, n_units = 16, n_treated = 5)
  d <- p$data
  d$id <- ifelse(d$unit == 2, 1 + .Machine$double.eps, d$unit)
  fit <- suppressWarnings(trop(d, "y", "d", "id", "time", lambda = l_r, se = "none",
                               control = small_control()))
  k <- unit_keys(fit$units, length(fit$units))
  w <- seq_along(k)
  a <- trop_effects(fit, by = "overall", unit_weights = setNames(w, k), n_boot = 0)
  b <- trop_effects(fit, by = "overall", unit_weights = w, n_boot = 0)
  expect_identical(a$estimate, b$estimate)
  expect_error(trop_effects(fit, unit_weights = setNames(w, as.character(fit$units))),
               "must be different")
})

test_that("the weights of units outside of the groups do not stop the call", {
  p <- small_panel(seed = 35, n_units = 16, n_treated = 5)
  fit <- trop_matrix(p$Y, p$W, lambda = l_r, se = "none", control = small_control())
  tr <- which(rowSums(p$W) > 0)
  grp <- matrix(NA_integer_, nrow(p$W), ncol(p$W))
  grp[tr[1], ] <- 1L
  w <- rep(1, nrow(p$W))
  w[tr[1]] <- 1e-170
  w[tr[2]] <- 1e50
  e <- trop_effects(fit, by = grp, unit_weights = w, n_boot = 0)
  expect_equal(e$estimate, mean(fit$tau[tr[1], ][is.finite(fit$tau[tr[1], ])]))
  expect_equal(trop_cells(fit, by = grp, unit_weights = w)$share, 1)
  expect_equal(nrow(trop_cells(fit, window = c(50, 60), unit_weights = w)), 0L)
})

test_that("a cell without a row between two runs with different starts does not join them", {
  set.seed(71)
  N <- 8
  TT <- 8
  ids <- sprintf("u%d", seq_len(N))
  d <- expand.grid(unit = ids, time = seq_len(TT), stringsAsFactors = FALSE)
  i1 <- d$unit == "u1"
  i2 <- d$unit == "u2"
  d$treated <- 0L
  d$treated[i1] <- as.integer(d$time[i1] %in% c(2, 3, 5:8))
  d$treated[i2] <- as.integer(d$time[i2] >= 5)
  d$event_time <- NA_real_
  d$event_time[i1] <- c(NA, 0, 1, NA, 0:3)
  d$event_time[i2] <- d$time[i2] - 5
  d$y <- rnorm(N)[match(d$unit, ids)] + rnorm(TT)[d$time] + 2 * d$treated
  d$y[i1 & d$time == 4] <- NA
  fit_d <- function(x) {
    suppressWarnings(trop(x, "y", "treated", "unit", "time", event_time = "event_time",
                          lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                          control = small_control()))
  }
  full <- fit_d(d)
  sparse <- fit_d(d[!(i1 & d$time == 4), ])
  expect_identical(sparse$tau, full$tau)
  expect_identical(sparse$panel$W, full$panel$W)
  co <- trop_cells(sparse, by = "cohort")
  expect_equal(co$cohort, c(2L, 5L))
  expect_equal(co$n_cells, c(2L, 8L))
})

test_that("the global cross-validation uses the drawn placebo sets of several units", {
  p <- small_panel(seed = 81, n_units = 16, n_treated = 5)
  g <- trop_grid(time = 0, unit = 0.5, nn = Inf)
  args <- list(p$Y, p$W, method = "global", grid = g, n_placebo = 6L, placebo_size = 3L,
               seed = 9, se = "none", control = tight_control())
  fit <- suppressWarnings(do.call(trop_matrix, args))
  plan <- do.call(trop_prepare, args)
  expect_true(all(apply(plan$sets, 1, anyDuplicated) == 0L))
  Yc <- plan$Yc
  TT <- ncol(Yc)
  k <- plan$design$k
  ate <- apply(plan$sets + 1L, 1, function(j) {
    Wp <- matrix(0L, nrow(Yc), TT)
    Wp[j, (TT - k + 1L):TT] <- 1L
    suppressWarnings(trop_matrix(Yc, Wp, method = "global", lambda = c(time = 0, unit = 0.5, nn = Inf),
                                 se = "none", control = tight_control()))$estimate
  })
  expect_equal(fit$cv$score, sqrt(mean(ate^2)), tolerance = 1e-8)
  expect_equal(fit$cv$n_cells, 6L * 3L * k)
  expect_equal(fit$cv$n_valid_chosen, fit$cv$n_cells)
})

test_that("a placebo unit without held-out outcomes does not change the weights", {
  N <- 13L
  TT <- 7L
  W <- matrix(0L, N, TT)
  W[1:5, 5:7] <- 1L
  Y <- outer(c(rep(0, 5), 0, 60, 0, 0, 1, 1, 1, 1), (0:6) / 2, "+")
  Y[c(8, 9), 5:7] <- Y[c(8, 9), 5:7] + 4
  Y[7, 5:7] <- NA_real_
  Y[1:5, 5:7] <- Y[1:5, 5:7] + 2
  g <- trop_grid(time = 0, unit = 1, nn = Inf)
  args <- list(Y, W, method = "global", grid = g, n_placebo = 1L, placebo_size = 2L,
               se = "none", control = tight_control())
  # Find a seed that draws the control units 6 and 7. Unit 7 has no outcome in
  # the held-out periods.
  seed <- NA
  for (s in 1:200) {
    plan <- do.call(trop_prepare, c(args, seed = s))
    if (setequal(c(plan$sets) + 1L, 1:2)) {
      seed <- s
      break
    }
  }
  expect_false(is.na(seed))
  fit <- suppressWarnings(do.call(trop_matrix, c(args, seed = seed)))
  Wp <- matrix(0L, nrow(plan$Yc), TT)
  Wp[1:2, 5:7] <- 1L
  ref <- suppressWarnings(trop_matrix(plan$Yc, Wp, method = "global",
                                      lambda = c(time = 0, unit = 1, nn = Inf), se = "none",
                                      control = tight_control()))
  expect_equal(ref$n_treated_units, 1L)
  expect_equal(fit$cv$score, abs(ref$estimate), tolerance = 1e-10)
  expect_equal(fit$cv$trace$n_valid, 3L)
})

test_that("the means of very large effects are finite, and the checkpoint reloads", {
  N <- 6L
  TT <- 7L
  W <- matrix(0L, N, TT)
  W[1:2, 4:7] <- 1L
  Y <- 6e307 * W
  fit <- trop_matrix(Y, W, lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                     control = small_control())
  expect_equal(fit$estimate, 6e307)
  ck <- tempfile(fileext = ".rds")
  on.exit(unlink(ck), add = TRUE)
  run <- function() {
    trop_effects(fit, by = list(event = "event_time", all = "overall"),
                 boot_rows = rbind(seq_len(N), seq_len(N)), checkpoint = ck)
  }
  a <- run()
  expect_equal(a$event$estimate, rep(6e307, 4))
  expect_equal(a$all$estimate, 6e307)
  expect_equal(a$all$n_draws, 2L)
  expect_true(all(is.finite(attr(a$event, "draws"))))
  b <- run()
  expect_identical(attr(b$all, "draws"), attr(a$all, "draws"))
})

# Round 4 of the review.

test_that("event clocks decide a cell without a row only when a new run starts after it", {
  N <- 8L
  TT <- 9L
  ids <- sprintf("u%d", seq_len(N))
  z <- function(d, t) d$unit == "u1" & d$time %in% t
  d <- expand.grid(unit = ids, time = seq_len(TT), stringsAsFactors = FALSE)
  d$w <- as.integer(z(d, c(2:4, 6:9)))
  d$event <- ifelse(d$w == 1L, ifelse(d$time <= 4L, d$time - 2L, d$time - 6L), NA_integer_)
  d$y <- match(d$unit, ids) + d$time / 4 + ifelse(d$w == 1L, 1 + d$event, 0)
  d$y[z(d, 5L)] <- NA_real_
  ff <- function(x) {
    suppressWarnings(trop(x, "y", "w", "unit", "time", event_time = "event",
                          lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                          control = small_control()))
  }
  # The treated row just before the gap has no outcome and no event time. The
  # rows after the gap start a new run in period 6.
  d$y[z(d, 4L)] <- NA_real_
  d$event[z(d, 4L)] <- NA_integer_
  ref <- ff(d)
  sparse <- ff(d[!z(d, 5L), ])
  expect_identical(sparse$panel$W, ref$panel$W)
  expect_identical(sparse$tau, ref$tau)
  expect_equal(unname(sparse$panel$W[1, 5]), 0L)

  # A clock that puts the start before the gap must agree with the rows before
  # the gap. Here the treatment starts in period 3, but the clock after the gap
  # says period 2.
  b <- expand.grid(unit = ids, time = seq_len(TT), stringsAsFactors = FALSE)
  b$w <- as.integer(z(b, 3:9))
  b$event <- ifelse(b$w == 1L, b$time - 3L, NA_integer_)
  b$y <- match(b$unit, ids) + b$time / 4 + ifelse(b$w == 1L, 1 + b$event, 0)
  wrong <- b[!z(b, 5L), ]
  wrong$event[z(wrong, 6:9)] <- 4:7
  expect_error(ff(wrong), "do not increase by 1")
})

test_that("the estimate of a group does not depend on the weights of other groups", {
  N <- 9L
  TT <- 8L
  W <- matrix(0L, N, TT)
  W[1:5, 5:8] <- 1L
  Y <- matrix(0, N, TT)
  Y[1, 5:8] <- 2e-140
  Y[2:5, 5:8] <- 1
  fit <- trop_matrix(Y, W, lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                     control = small_control())
  a <- b <- matrix(NA_integer_, N, TT)
  a[1, 5:8] <- 1L
  b[2, 5:8] <- 1L
  ident <- rbind(seq_len(N), seq_len(N))
  for (tiny in c(1e-190, 1e-210)) {
    uw <- c(tiny, rep(1, N - 1L))
    for (rule in c("cell", "unit")) {
      one <- trop_effects(fit, by = a, weights = rule, unit_weights = uw, boot_rows = ident)
      two <- trop_effects(fit, by = list(a = a, b = b), weights = rule, unit_weights = uw,
                          boot_rows = ident)
      expect_equal(one$estimate, 2e-140)
      expect_identical(two$a$estimate, one$estimate)
      expect_identical(unname(attr(two$a, "draws")), unname(attr(one, "draws")))
      expect_equal(unname(attr(one, "draws"))[, 1], rep(2e-140, 2))
    }
  }
  # In one group, the range of the weights is checked.
  uw <- c(1e-210, rep(1, N - 1L))
  expect_error(trop_effects(fit, by = "overall", unit_weights = uw, n_boot = 0), "too different")
})

test_that("unit names do not depend on the decimal mark of the session", {
  op <- options(OutDec = ",")
  on.exit(options(op), add = TRUE)
  expect_equal(unit_keys(c(1.125, 2.5, 3), 3), c("1.125", "2.5", "3"))
  p <- small_panel(seed = 37, n_units = 14, n_treated = 5)
  d <- p$data
  d$id <- d$unit + 0.125
  rows <- trop_boot_rows(p$W, 8, seed = 2)
  fit <- suppressWarnings(trop(d, "y", "d", "id", "time", lambda = l_r, boot_rows = rows,
                               control = small_control()))
  ct <- trop_contrast(list(a = fit, b = fit), c(a = 1, b = -1))
  expect_equal(ct$std.error, 0)
  w <- setNames(seq_along(fit$units), unit_keys(fit$units, length(fit$units)))
  expect_s3_class(trop_effects(fit, by = "overall", unit_weights = w, n_boot = 0), "trop_effects")
})

test_that("the comparison of the draws does not depend on the order of the objects", {
  p <- small_panel(seed = 39, n_units = 14, n_treated = 5)
  n <- nrow(p$Y)
  id <- seq_len(n)
  rows <- trop_boot_rows(p$W, 10, seed = 6)
  pm <- c(2L, 1L, id[-(1:2)])
  mp <- matrix(match(rows, pm), nrow(rows))
  mf <- function(names, idx = id, r = rows) {
    y <- p$Y[idx, , drop = FALSE]
    rownames(y) <- names
    trop_matrix(y, p$W[idx, , drop = FALSE], lambda = l_r, boot_rows = r,
                control = small_control())
  }
  A <- mf(as.character(id))
  B <- mf(as.character(id)[pm], pm, mp)
  U <- mf(NULL, pm, mp)
  # A and B have the same units in each draw, with other row numbers. U has no
  # names. Thus, U must have the same row numbers as A and B, which is false.
  expect_error(suppressWarnings(trop_stepdown(list(A = A, B = B, U = U))), "same units")
  expect_error(suppressWarnings(trop_stepdown(list(U = U, B = B, A = A))), "same units")
  expect_equal(nrow(trop_stepdown(list(A = A, B = B))), 2L)
  expect_warning(st <- trop_stepdown(list(B = B, U = U)), "no valid unit names: U")
  expect_equal(nrow(st), 2L)
})
