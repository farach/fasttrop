# Council condition C2: tuning and inference correctness.

test_that("leave-one-out CV holds each cell out of its own fit and its own distances", {
  p <- small_panel(seed = 13, n_units = 10, n_periods = 8, n_treated = 3, n_post = 2)
  inp <- fasttrop:::as_panel_inputs(p$Y, p$W)
  l <- c(time = 0.5, unit = 0.8, nn = 0.5)
  ctl <- tight_control()
  pool <- which(inp$W == 0L) - 1L
  cells <- as.integer(sort(pool[c(3, 17, 29, 40, 52)]))
  kernel <- .Call(fasttrop:::C_local_loocv, inp$Y, inp$W, inp$obs, l, cells, unclass(ctl))
  # Reference: mark each held-out cell as treated and estimate its effect.
  manual <- vapply(cells, function(c0) {
    i <- c0 %% nrow(inp$Y) + 1L
    t <- c0 %/% nrow(inp$Y) + 1L
    Wm <- inp$W
    Wm[i, t] <- 1L
    fit <- suppressWarnings(trop_matrix(inp$Y, Wm, lambda = l, se = "none", control = ctl))
    fit$tau[i, t]
  }, double(1))
  expect_equal(kernel$n_valid, length(cells))
  expect_equal(kernel$Q, sum(manual^2), tolerance = 1e-8)
})

test_that("placebo CV is the RMSE of placebo ATEs on untreated units", {
  p <- small_panel(seed = 17, n_units = 14, n_periods = 9, n_treated = 3, n_post = 3)
  never <- which(rowSums(p$W) == 0)
  Yc <- p$Y[never, ]
  k <- 3L
  sets <- matrix(c(1L, 4L, 6L, 2L), ncol = 1L) - 1L
  l <- c(time = 0.2, unit = 0.6, nn = 0.5)
  ctl <- tight_control()
  oc <- matrix(1L, nrow(Yc), ncol(Yc))
  Wz <- matrix(0L, nrow(Yc), ncol(Yc))
  kernel <- .Call(fasttrop:::C_global_placebo, Yc, oc, Wz, k, sets, l, unclass(ctl))
  ates <- vapply(sets[, 1] + 1L, function(u) {
    Wp <- Wz
    Wp[u, (ncol(Yc) - k + 1L):ncol(Yc)] <- 1L
    trop_matrix(Yc, Wp, method = "global", lambda = l, se = "none", control = ctl)$estimate
  }, double(1))
  expect_equal(kernel$rmse, sqrt(mean(ates^2)), tolerance = 1e-8)
})

test_that("the local search reports the score of the chosen grid point", {
  p <- small_panel(seed = 19, n_units = 10, n_periods = 8, n_post = 2)
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 1), nn = c(0.5, Inf))
  fit <- trop_matrix(p$Y, p$W, grid = g, cv_cells = 20, se = "none", seed = 2,
                     control = small_control())
  tr <- fit$cv$trace
  key <- paste(tr$time, tr$unit, tr$nn)
  b <- fit$lambda
  expect_equal(fit$cv$score, tr$score[key == paste(b[["time"]], b[["unit"]], b[["nn"]])])
  expect_true(all(tr$time %in% g$time) && all(tr$unit %in% g$unit) && all(tr$nn %in% g$nn))
  expect_equal(fit$cv$n_evals, nrow(tr))
  expect_equal(fit$cv$n_cells, 20L)
})

test_that("the global search reaches a fixed point of the coordinate cycle", {
  p <- small_panel(seed = 21, n_units = 16, n_periods = 9, n_post = 3)
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 1), nn = c(0.5, Inf))
  fit <- trop_matrix(p$Y, p$W, method = "global", grid = g, n_placebo = 20, se = "none",
                     seed = 4, control = small_control())
  expect_true(fit$cv$fixed_point)
  expect_equal(fit$cv$n_sets, 20L)
  tr <- fit$cv$trace
  b <- fit$lambda
  # At a fixed point no single coordinate change on the grid improves the score.
  axis_scores <- function(axis, v) {
    pt <- b
    pt[[axis]] <- v
    tr$score[tr$time == pt[["time"]] & tr$unit == pt[["unit"]] & tr$nn == pt[["nn"]]]
  }
  for (axis in c("time", "unit", "nn")) {
    for (v in g[[axis]]) {
      s <- axis_scores(axis, v)
      expect_length(s, 1L)
      expect_gte(s, fit$cv$score)
    }
  }
})

test_that("bootstrap draws equal re-estimation on the resampled panel", {
  p <- small_panel(seed = 23)
  l <- c(time = 0.3, unit = 0.5, nn = 0.5)
  for (method in c("local", "global")) {
    plan <- fasttrop:::trop_prepare(p$Y, p$W, method = method, lambda = l, n_boot = 4,
                                    seed = 9, control = tight_control())
    fit <- fasttrop:::trop_run(plan)
    for (b in 1:4) {
      rows <- plan$rows[b, ] + 1L
      direct <- suppressWarnings(trop_matrix(p$Y[rows, ], p$W[rows, ], method = method,
                                             lambda = l, se = "none", control = tight_control()))
      expect_equal(fit$boot$draws[b], direct$estimate, tolerance = 1e-8,
                   label = paste(method, "draw", b))
    }
  }
})

test_that("resampling follows Algorithm 3, with a controls-only option", {
  p <- small_panel(seed = 29)
  tr <- which(rowSums(p$W) > 0)
  never <- setdiff(seq_len(nrow(p$W)), tr)
  l <- c(time = 0, unit = 0, nn = Inf)
  prep <- function(resample) {
    fasttrop:::trop_prepare(p$Y, p$W, lambda = l, n_boot = 30, seed = 1, resample = resample,
                            control = small_control())$rows + 1L
  }
  rs <- prep("stratified")
  rc <- prep("controls")
  slots <- seq_along(tr)
  expect_true(all(rs[, slots] %in% tr))
  expect_true(all(rs[, -slots] %in% never))
  expect_true(any(apply(rs[, slots, drop = FALSE], 1, anyDuplicated) > 0))
  expect_true(all(apply(rc[, slots, drop = FALSE], 1, identical, tr)))
  expect_true(all(rc[, -slots] %in% never))
})

test_that("a seed makes fits reproducible and leaves the global RNG untouched", {
  p <- small_panel(seed = 31)
  l <- c(time = 0.2, unit = 0.4, nn = 1)
  set.seed(100)
  before <- .Random.seed
  f1 <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 15, seed = 42, control = small_control())
  expect_identical(.Random.seed, before)
  f2 <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 15, seed = 42, control = small_control())
  expect_identical(f1$boot$draws, f2$boot$draws)
  expect_equal(f1$std.error, stats::sd(f1$boot$draws))
  expect_identical(f1$ci_type, "normal")
  expect_equal(unname(confint(f1)[1, ]),
               f1$estimate + c(-1, 1) * stats::qnorm(0.975) * f1$std.error)
  f3 <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 15, seed = 42, ci = "percentile",
                    control = small_control())
  expect_identical(f3$boot$draws, f1$boot$draws)
  expect_equal(unname(confint(f3)[1, ]),
               unname(stats::quantile(f3$boot$draws, c(0.025, 0.975))))
})
