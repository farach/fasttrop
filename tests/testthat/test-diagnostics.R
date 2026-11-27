# Tests for the diagnostics of the fits, the SVD failure path, and the
# distance rows of the C code.

l_fin <- c(time = 0.3, unit = 0.5, nn = 1)

# The local estimate on the resampled panel of one draw.
draw_estimate <- function(Y, W, rows_b, lambda, ctl) {
  Yb <- Y[rows_b, , drop = FALSE]
  Wb <- W[rows_b, , drop = FALSE] * 1L
  inp <- fasttrop:::as_panel_inputs(Yb, Wb)
  .Call(fasttrop:::C_local_estimate, inp$Y, inp$W, inp$obs, NULL, lambda, unclass(ctl))
}

test_that("the counts of each draw equal the counts of a fit on the resampled panel", {
  p <- small_panel(seed = 101, n_units = 16, n_treated = 5)
  ctl <- trop_control(threads = 2L, max_iter = 1L, max_inner = 1L)
  rows <- trop_boot_rows(p$W, 6, seed = 3)
  w <- collect_warnings(fit <- trop_matrix(p$Y, p$W, lambda = l_fin, boot_rows = rows, control = ctl))
  d <- fit$boot$diagnostics
  expect_named(d, c("draw", "n_used", "n_skipped", "n_nonconverged", "n_failed", "n_dropped"))
  expect_equal(nrow(d), 6L)
  for (b in seq_len(nrow(rows))) {
    ref <- draw_estimate(p$Y, p$W, rows[b, ], l_fin, ctl)
    expect_equal(d$n_used[b], ref$n_used)
    expect_equal(d$n_nonconverged[b], ref$n_nonconverged)
    expect_equal(fit$boot$draws[b], ref$att, tolerance = 1e-12)
  }
  expect_true(all(d$n_nonconverged > 0L))
  expect_true(any(grepl("bootstrap draws, [0-9]+ fit\\(s\\) did not converge", w)))
  expect_output(print(summary(fit)), "Bootstrap fits: [0-9]+ fit\\(s\\) in 6 draw\\(s\\) did not converge")
})

test_that("a draw with duplicated units uses the distances of the original units", {
  p <- small_panel(seed = 102, n_units = 14, n_treated = 5)
  ctl <- tight_control()
  rows <- rbind(c(1, 1, 2, 3, 3, 3, 7:14), c(5, 4, 4, 1, 2, 2, 6, 6, 6, 9:13))
  l <- c(time = 0.4, unit = 1.2, nn = 0.7)
  fit <- trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows, control = ctl)
  for (b in 1:2) {
    ref <- draw_estimate(p$Y, p$W, rows[b, ], l, ctl)
    expect_equal(fit$boot$draws[b], ref$att, tolerance = 1e-12)
  }
})

test_that("the second SVD routine gives the same result and a double failure fails the fit", {
  p <- small_panel(seed = 103, n_units = 16, n_treated = 5)
  rows <- trop_boot_rows(p$W, 6, seed = 2)
  ok <- trop_matrix(p$Y, p$W, lambda = l_fin, boot_rows = rows, control = small_control())
  with_fault <- function(code) {
    ctl <- small_control()
    ctl$fault <- code
    ctl
  }
  f1 <- trop_matrix(p$Y, p$W, lambda = l_fin, boot_rows = rows, control = with_fault(1L))
  expect_equal(f1$estimate, ok$estimate, tolerance = 1e-10)
  expect_equal(f1$boot$draws, ok$boot$draws, tolerance = 1e-10)
  expect_equal(f1$n_failed, 0L)

  w <- collect_warnings(f3 <- trop_matrix(p$Y, p$W, lambda = l_fin, boot_rows = rows,
                                          control = with_fault(3L)))
  expect_true(is.na(f3$estimate))
  expect_true(all(is.na(f3$tau)))
  expect_equal(f3$n_failed, sum(p$W))
  expect_equal(f3$n_nonconverged, 0L)
  expect_true(all(f3$boot$diagnostics$n_used == 0L))
  expect_true(all(f3$boot$diagnostics$n_failed > 0L))
  expect_true(any(grepl("^[0-9]+ fit\\(s\\) failed", w)))
  expect_true(any(grepl("bootstrap draws, [0-9]+ fit\\(s\\) failed", w)))
  expect_output(print(summary(f3)), "Failed fits: 15")

  # With nn = Inf, the fit has no SVD.
  inf <- trop_matrix(p$Y, p$W, lambda = c(time = 0.3, unit = 0.5, nn = Inf), se = "none",
                     control = with_fault(3L))
  expect_equal(inf$n_failed, 0L)
  expect_true(is.finite(inf$estimate))

  # Bit 4: the faults occur only in the draws with an odd index in the C code.
  w7 <- collect_warnings(f7 <- trop_matrix(p$Y, p$W, lambda = l_fin, boot_rows = rows,
                                           control = with_fault(7L)))
  expect_identical(f7$estimate, ok$estimate)
  expect_true(all(is.na(f7$boot$draws[c(2, 4, 6)])))
  expect_identical(f7$boot$draws[c(1, 3, 5)], ok$boot$draws[c(1, 3, 5)])
  expect_equal(f7$boot$n_failed, 3L)
  e7 <- suppressWarnings(trop_effects(f7, by = "event_time"))
  expect_true(all(e7$n_draws == 3L))
  expect_equal(attr(e7, "n_complete"), 3L)

  gl <- suppressWarnings(trop_matrix(p$Y, p$W, method = "global", lambda = l_fin, se = "none",
                                     control = with_fault(3L)))
  expect_equal(gl$n_failed, 1L)
  expect_equal(gl$n_nonconverged, 0L)
  expect_true(is.na(gl$estimate))
  expect_error(fasttrop:::trop_fit_weighted(p$Y, matrix(1, nrow(p$Y), ncol(p$Y)), 1,
                                            control = with_fault(3L)),
               "singular value decomposition failed")
})

test_that("the cross-validation counts the fits of each grid point", {
  p <- small_panel(seed = 105)
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 1), nn = c(1, Inf))
  ctl <- trop_control(threads = 2L, max_iter = 2L, max_inner = 1L)
  w <- collect_warnings(fit <- trop_matrix(p$Y, p$W, grid = g, cv_cells = 30, se = "none",
                                           seed = 4, control = ctl))
  tr <- fit$cv$trace
  expect_named(tr, c("time", "unit", "nn", "score", "n_valid", "n_nonconverged", "n_failed"))
  expect_equal(fit$cv$n_nonconverged, sum(tr$n_nonconverged))
  expect_true(all(tr$n_nonconverged[is.infinite(tr$nn)] == 0L))
  expect_true(any(tr$n_nonconverged > 0L))
  expect_true(any(grepl("In the cross-validation, [0-9]+ leave-one-out fit\\(s\\) did not converge", w)))
  plan <- fasttrop:::trop_prepare(p$Y, p$W, grid = g, cv_cells = 30, se = "none", seed = 4,
                                  control = ctl)
  for (k in seq_len(nrow(tr))) {
    cv <- .Call(fasttrop:::C_local_loocv, plan$p$Y, plan$p$W, plan$p$obs, NULL,
                unname(unlist(tr[k, 1:3])), plan$cells, unclass(ctl))
    expect_equal(tr$score[k], cv$Q, tolerance = 0)
    expect_equal(tr$n_valid[k], cv$n_valid)
    expect_equal(tr$n_nonconverged[k], cv$n_nonconverged)
  }
  expect_output(print(summary(fit)), "Tuning fits: [0-9]+ did not converge, 0 failed")
})

# A panel in which only control unit 12 has an outcome in the last period.
# A draw without unit 12 cannot identify the treated cells of that period.
sparse_panel <- function() {
  set.seed(11)
  N <- 12
  TT <- 6
  W <- matrix(0L, N, TT)
  W[1:4, 5:6] <- 1L
  Y <- outer(rnorm(N), rep(1, TT)) + outer(rep(1, N), rnorm(TT)) + W +
    matrix(rnorm(N * TT, 0, 0.1), N, TT)
  Y[5:11, 6] <- NA
  list(Y = Y, W = W)
}

test_that("n_dropped counts the cells that a draw cannot identify, and the band uses complete draws", {
  s <- sparse_panel()
  rows <- trop_boot_rows(s$W, 30, seed = 8)
  twfe <- c(time = 0, unit = 0, nn = Inf)
  w <- collect_warnings(fit <- trop_matrix(s$Y, s$W, lambda = twfe, boot_rows = rows,
                                           control = small_control()))
  expect_equal(fit$n_skipped, 0L)
  has12 <- apply(rows, 1L, function(r) 12L %in% r)
  expect_true(any(!has12) && any(has12))
  d <- fit$boot$diagnostics
  expect_equal(d$n_dropped, ifelse(has12, 0L, 4L))
  expect_true(any(grepl("do not identify [0-9]+ treated cell", w)))

  ev <- suppressWarnings(trop_effects(fit, by = "event_time"))
  expect_equal(ev$n_draws, c(30L, sum(has12)))
  expect_equal(attr(ev, "n_complete"), sum(has12))
  # The band uses only the complete draws, for the centers, the standard
  # errors, and the critical value.
  D <- attr(ev, "draws")[has12, ]
  se_c <- apply(D, 2L, stats::sd)
  Z <- abs(sweep(sweep(D, 2L, colMeans(D)), 2L, se_c, "/"))
  crit <- unname(stats::quantile(apply(Z, 1L, max), 0.95, type = 7))
  expect_equal(attr(ev, "critical"), crit)
  expect_equal(ev$band.high - ev$estimate, crit * se_c, tolerance = 1e-12, ignore_attr = TRUE)
  se_all <- apply(attr(ev, "draws"), 2L, function(v) stats::sd(v[is.finite(v)]))
  expect_equal(ev$std.error, unname(se_all))
  expect_output(print(ev), sprintf("It uses the %d draws with a value for all rows", sum(has12)))
  expect_output(print(ev), "the band and the intervals use different draws")
  expect_warning(sd <- trop_stepdown(ev), "do not have a value for all estimates")
  expect_equal(attr(sd, "n_draws"), sum(has12))
  expect_output(print(sd), sprintf("Draws: %d of 30 have a value for all estimates", sum(has12)))
})

test_that("the fit records the C time of each step", {
  p <- small_panel(seed = 107)
  fit <- trop_matrix(p$Y, p$W, grid = trop_grid(time = c(0, 0.5), unit = 0, nn = Inf),
                     cv_cells = 20, n_boot = 5, seed = 1, control = small_control())
  tm <- fit$timing
  expect_equal(tm$kernel_seconds, tm$cv_seconds + tm$estimate_seconds + tm$bootstrap_seconds)
  expect_true(all(c(tm$cv_seconds, tm$estimate_seconds, tm$bootstrap_seconds) > 0))
  expect_output(print(summary(fit)), "C time by step: tuning")
})
