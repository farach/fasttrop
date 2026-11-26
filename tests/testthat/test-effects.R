# Tests for trop_effects() and trop_stepdown().

twfe <- c(time = 0, unit = 0, nn = Inf)

# The first treated period of each unit. A never-treated unit has NA.
first_treated <- function(W) {
  s <- apply(W, 1L, function(w) if (any(w == 1L)) which(w == 1L)[1L] else NA_integer_)
  as.integer(s)
}

# The true mean effect of the treated cells in each group.
group_truth <- function(p, key) {
  k <- key[p$W == 1L]
  c(tapply(p$effect[p$W == 1L][!is.na(k)], k[!is.na(k)], mean))
}

# A panel with 2 covariates for the check of the draws.
cov_panel_for_effects <- function(seed = 5, N = 24, TT = 10) {
  set.seed(seed)
  a <- rnorm(N)
  b <- rnorm(TT)
  x1 <- matrix(rnorm(N * TT), N, TT) + a
  x2 <- matrix(runif(N * TT), N, TT) + outer(rep(1, N), b)
  W <- matrix(0L, N, TT)
  W[1:6, (TT - 2):TT] <- 1L
  Y <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + 2 * x1 - 1.5 * x2 + 0.7 * W +
    matrix(rnorm(N * TT, sd = 0.3), N, TT)
  X <- array(c(x1, x2), c(N, TT, 2), dimnames = list(NULL, NULL, c("x1", "x2")))
  list(Y = Y, W = W, X = X)
}

staggered_fit <- function(...) {
  p <- staggered_panel()
  fit <- suppressWarnings(trop_matrix(p$Y, p$W, lambda = twfe, se = "none",
                                      control = tight_control(), ...))
  list(p = p, fit = fit)
}

test_that("trop_effects recovers the true effects of each group", {
  s <- staggered_fit()
  p <- s$p
  start <- first_treated(p$W)
  event <- col(p$W) - start[row(p$W)]

  ev <- trop_effects(s$fit, by = "event_time")
  expect_s3_class(ev, "trop_effects")
  truth <- group_truth(p, event)
  expect_equal(ev$event_time, as.numeric(names(truth)))
  expect_equal(ev$estimate, unname(truth), tolerance = 1e-8)

  co <- trop_effects(s$fit, by = "cohort")
  expect_equal(co$cohort, c(5, 7, 9))
  expect_equal(co$estimate, unname(group_truth(p, matrix(start, nrow(p$W), ncol(p$W)))),
               tolerance = 1e-8)
  expect_equal(co$n_cells, c(3 * 8, 3 * 6, 3 * 4))

  tm <- trop_effects(s$fit, by = "time")
  expect_equal(tm$estimate, unname(group_truth(p, col(p$W) * 1L)), tolerance = 1e-8)

  un <- trop_effects(s$fit, by = "unit")
  expect_equal(un$unit, 1:9)
  expect_equal(un$estimate, unname(group_truth(p, row(p$W) * 1L)), tolerance = 1e-8)

  ov <- trop_effects(s$fit, by = "overall")
  expect_equal(nrow(ov), 1L)
  expect_equal(ov$estimate, s$fit$estimate, tolerance = 1e-10)
  expect_equal(ov$estimate, p$att, tolerance = 1e-8)
  expect_equal(ov$n_cells, sum(p$W))
  expect_false("std.error" %in% names(ov))
})

test_that("window limits the event times and a matrix by makes any groups", {
  s <- staggered_fit()
  ev <- trop_effects(s$fit, by = "event_time", window = c(0, 1))
  expect_equal(ev$event_time, c(0, 1))
  expect_equal(ev$n_cells, c(9, 9))
  expect_error(trop_effects(s$fit, by = "event_time", window = c(2, 1)), "window")
  expect_error(trop_effects(s$fit, by = "event_time", window = 1), "window")

  W <- s$p$W
  by <- matrix(NA_integer_, nrow(W), ncol(W))
  by[W == 1L] <- 1L
  by[1:3, ][W[1:3, ] == 1L] <- 2L
  by[, 1:6][W[, 1:6] == 1L] <- NA
  g <- trop_effects(s$fit, by = by)
  expect_identical(names(g)[1L], "group")
  expect_equal(g$n_cells, c(sum(!is.na(by) & by == 1L), sum(!is.na(by) & by == 2L)))
  expect_error(trop_effects(s$fit, by = by[-1, ]), "dimensions")
})

test_that("trop_effects uses the labels of a fit from a data frame", {
  p <- small_panel(seed = 71)
  d <- p$data
  d$unit <- paste0("u", d$unit)
  d$time <- 2000 + d$time
  fit <- trop(d, "y", "d", "unit", "time", lambda = c(time = 0.3, unit = 0.5, nn = 1),
              control = small_control())
  tm <- trop_effects(fit, by = "time")
  expect_true(all(tm$time %in% (2000 + 1:10)))
  un <- trop_effects(fit, by = "unit")
  expect_true(all(grepl("^u", un$unit)))
})

test_that("the draws of trop_effects equal the draws of the fit", {
  p <- cov_panel_for_effects()
  l <- c(time = 0.3, unit = 0.4, nn = 0.5)
  cases <- list(
    list(),
    list(X = p$X, covariate_method = "joint"),
    list(X = p$X, covariate_method = "residualize")
  )
  for (a in cases) {
    fit <- do.call(trop_matrix, c(list(p$Y, p$W, lambda = l, n_boot = 12, seed = 5,
                                       control = tight_control()), a))
    ov <- trop_effects(fit, by = "overall")
    expect_equal(unname(attr(ov, "draws")[, 1L]), fit$boot$draws, tolerance = 1e-8)
    expect_equal(ov$std.error, fit$std.error, tolerance = 1e-8)
    expect_equal(ov$estimate, fit$estimate, tolerance = 1e-10)
  }
})

test_that("trop_effects makes new draws without changing the random state", {
  p <- small_panel(seed = 73)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = small_control())
  expect_null(attr(trop_effects(fit, by = "overall"), "rows"))
  set.seed(9)
  before <- .Random.seed
  a <- trop_effects(fit, by = "time", n_boot = 10, seed = 5)
  expect_identical(.Random.seed, before)
  expect_equal(attr(a, "rows"), trop_boot_rows(p$W, 10, seed = 5))
  expect_true(all(c("std.error", "conf.low", "conf.high", "band.low", "band.high") %in% names(a)))
  expect_equal(dim(attr(a, "draws")), c(10L, nrow(a)))
  expect_true(is.finite(attr(a, "critical")))
  expect_equal(a$conf.high - a$estimate, a$estimate - a$conf.low, tolerance = 1e-10)

  b <- trop_effects(fit, by = "time", n_boot = 10, seed = 5, ci = "percentile")
  expect_equal(b$std.error, a$std.error)
  d <- attr(b, "draws")
  expect_true(all(b$conf.low >= apply(d, 2L, min) - 1e-12))
  expect_true(all(b$conf.high <= apply(d, 2L, max) + 1e-12))

  one <- trop_effects(fit, by = "overall", n_boot = 10, seed = 5)
  expect_true(is.na(one$band.low))
  expect_null(attr(trop_effects(fit, by = "overall", n_boot = 0), "draws"))
})

test_that("trop_effects uses boot_rows and strata, and rejects wrong input", {
  p <- small_panel(seed = 75)
  N <- nrow(p$W)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = small_control())
  rows <- trop_boot_rows(p$W, 8, seed = 2)
  a <- trop_effects(fit, by = "overall", boot_rows = rows)
  expect_equal(attr(a, "rows"), rows)
  strata <- rep(c("a", "b"), length.out = N)
  b <- trop_effects(fit, by = "overall", n_boot = 8, strata = strata, seed = 2)
  expect_equal(attr(b, "rows"), trop_boot_rows(p$W, 8, strata = strata, seed = 2))
  expect_error(trop_effects(fit, by = "overall", boot_rows = rows, strata = strata), "not both")
  expect_error(trop_effects(fit, by = "overall", conf_level = 1), "conf_level")
  expect_error(trop_effects(list()), "trop")
  nopanel <- fit
  nopanel$panel <- NULL
  expect_error(trop_effects(nopanel), "does not contain the panel")
  gl <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none", control = small_control())
  expect_error(trop_effects(gl), "method = \"local\"")
})

test_that("trop_effects respects the target of the fit", {
  p <- small_panel(seed = 77)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  tgt <- matrix(TRUE, nrow(p$W), ncol(p$W))
  tgt[, ncol(p$W)] <- FALSE
  fit <- trop_matrix(p$Y, p$W, lambda = l, se = "none", target = tgt, control = small_control())
  ov <- trop_effects(fit, by = "overall")
  expect_equal(ov$estimate, fit$estimate, tolerance = 1e-10)
  expect_equal(ov$n_cells, fit$n_treated_cells)
})

test_that("the print method shows the band and the as.data.frame method drops the attributes", {
  p <- small_panel(seed = 79)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 10, seed = 1, control = small_control())
  ev <- trop_effects(fit, by = "event_time")
  expect_output(print(ev), "TROP effects by event time")
  expect_output(print(ev), "sup-t")
  df <- as.data.frame(ev)
  expect_identical(class(df), "data.frame")
  expect_null(attr(df, "draws"))
  expect_equal(nrow(df), nrow(ev))
})

test_that("trop_stepdown gives adjusted p-values that are ordered", {
  p <- small_panel(seed = 81, n_units = 20, n_treated = 6)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  rows <- trop_boot_rows(p$W, 60, seed = 3)
  fit <- trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows, control = small_control())
  ev <- trop_effects(fit, by = "event_time")
  sd <- trop_stepdown(ev)
  expect_named(sd, c("name", "estimate", "std.error", "t", "p_raw", "p_adj"))
  expect_equal(nrow(sd), nrow(ev))
  expect_true(all(sd$p_adj >= sd$p_raw - 1e-12))
  expect_true(all(sd$p_adj > 0 & sd$p_adj <= 1))
  ord <- order(abs(sd$t), decreasing = TRUE)
  expect_false(is.unsorted(sd$p_adj[ord]))

  # A named list of fits with the same rows. A shift of the outcome does not change the estimate.
  fit2 <- trop_matrix(p$Y + 1, p$W, lambda = l, boot_rows = rows, control = small_control())
  both <- trop_stepdown(list(a = fit, b = fit2))
  expect_equal(both$name, c("a", "b"))
  expect_equal(both$estimate[2L], both$estimate[1L], tolerance = 1e-6)
  mixed <- trop_stepdown(list(fit = fit, ev = ev))
  expect_equal(nrow(mixed), 1L + nrow(ev))
})

test_that("trop_stepdown rejects objects with different rows or no bootstrap", {
  p <- small_panel(seed = 83)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  a <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 8, seed = 1, control = small_control())
  b <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 8, seed = 2, control = small_control())
  none <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = small_control())
  expect_error(trop_stepdown(list(a, b)), "same bootstrap rows")
  expect_error(trop_stepdown(none), "must have a bootstrap")
  expect_error(trop_stepdown(1), "`x` must be")
  expect_error(trop_stepdown(list(a, 1)), "Each element")
})

test_that("anticipation gives event time 0 to the first period after the anticipation", {
  p <- small_panel(seed = 91, n_units = 20, n_treated = 6, design = "staggered")
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  start <- first_treated(p$W)
  W <- p$W
  tgt <- matrix(TRUE, nrow(W), ncol(W))
  for (i in which(!is.na(start) & start > 1L)) {
    W[i, start[i] - 1L] <- 1L
    tgt[i, start[i] - 1L] <- FALSE
  }
  fit <- suppressWarnings(trop_matrix(p$Y, W, lambda = l, se = "none", target = tgt,
                                      control = small_control()))
  ev <- trop_effects(fit, by = "event_time", anticipation = 1)
  ok <- !is.na(start) & start > 1L
  expect_equal(min(ev$event_time), 0)
  expect_equal(ev$n_cells[ev$event_time == 0], sum(ok))
  co <- trop_effects(fit, by = "cohort", anticipation = 1)
  expect_setequal(co$cohort, unique(start[ok]))
  expect_error(trop_effects(fit, anticipation = -1), "anticipation")
  expect_error(trop_effects(fit, anticipation = 1.5), "anticipation")
})

test_that("treated cells that the target excludes make no empty group", {
  p <- small_panel(seed = 93, n_units = 30, n_treated = 8, design = "staggered")
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  start <- first_treated(p$W)
  tgt <- matrix(TRUE, nrow(p$W), ncol(p$W))
  for (i in which(!is.na(start))) tgt[i, start[i]] <- FALSE
  rows <- trop_boot_rows(p$W, 20, seed = 4)
  fit <- suppressWarnings(trop_matrix(p$Y, p$W, lambda = l, target = tgt, boot_rows = rows,
                                      control = small_control()))
  ev <- trop_effects(fit, by = "event_time")
  expect_false(0 %in% ev$event_time)
  expect_true(all(ev$n_cells > 0L))
  expect_true(all(is.finite(ev$estimate)))
  expect_equal(nrow(trop_stepdown(ev)), nrow(ev))
})

test_that("a list by gives the results of separate calls from one set of draws", {
  p <- small_panel(seed = 95, n_units = 20, n_treated = 6, design = "staggered")
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  rows <- trop_boot_rows(p$W, 15, seed = 6)
  fit <- suppressWarnings(trop_matrix(p$Y, p$W, lambda = l, se = "none", control = small_control()))
  m <- matrix(rep(1:2, length.out = length(p$W)), nrow(p$W))
  both <- trop_effects(fit, by = list(event = "event_time", "cohort", half = m), boot_rows = rows)
  expect_named(both, c("event", "cohort", "half"))
  ev <- trop_effects(fit, by = "event_time", boot_rows = rows)
  co <- trop_effects(fit, by = "cohort", boot_rows = rows)
  expect_equal(both$event, ev, tolerance = 1e-12)
  expect_equal(both$cohort, co, tolerance = 1e-12)
  expect_s3_class(both$half, "trop_effects")
  expect_equal(nrow(trop_stepdown(both)), nrow(ev) + nrow(co) + nrow(both$half))
  expect_error(trop_effects(fit, by = list()), "at least 1")
  expect_error(trop_effects(fit, by = list(m)), "Give a name")
  expect_error(trop_effects(fit, by = list("cohort", cohort = "time")), "must be different")
  expect_error(trop_effects(fit, by = list(1)), "one group name or a matrix")
})
