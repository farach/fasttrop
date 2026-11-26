# Tests for the target cells, the strata, the bootstrap rows, the demeaning, and
# the scaled grid.

last_period_off <- function(W) {
  tgt <- matrix(TRUE, nrow(W), ncol(W))
  tgt[, ncol(W)] <- FALSE
  tgt
}

test_that("a target gives the same fit as NA outcomes in the cells outside of it", {
  p <- small_panel(seed = 41)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  tgt <- last_period_off(p$W)
  fit <- trop_matrix(p$Y, p$W, lambda = l, se = "none", target = tgt, control = tight_control())
  Ym <- p$Y
  Ym[p$W == 1L & !tgt] <- NA
  manual <- trop_matrix(Ym, p$W, lambda = l, se = "none", control = tight_control())
  expect_equal(fit$estimate, manual$estimate, tolerance = 1e-10)
  expect_equal(fit$tau, manual$tau, tolerance = 1e-10)
  expect_equal(fit$n_excluded_cells, sum(p$W == 1L & !tgt))
  expect_equal(fit$n_treated_cells, sum(p$W == 1L & tgt))
  all_cells <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = tight_control())
  expect_equal(all_cells$n_excluded_cells, 0L)
  expect_false(isTRUE(all.equal(fit$estimate, all_cells$estimate)))
})

test_that("a target works with a bootstrap and with a staggered design", {
  p <- staggered_panel()
  tgt <- matrix(TRUE, nrow(p$W), ncol(p$W))
  tgt[, 10:12] <- FALSE
  l <- c(time = 0, unit = 0, nn = Inf)
  fit <- suppressWarnings(trop_matrix(p$Y, p$W, lambda = l, target = tgt, n_boot = 20, seed = 3,
                                      control = small_control()))
  ok <- p$W == 1L & tgt
  # The data have no noise and parallel trends. Thus, TWFE imputation is exact.
  expect_equal(fit$estimate, mean(p$effect[ok]), tolerance = 1e-8)
  expect_equal(fit$n_treated_cells, sum(ok))
  expect_true(is.finite(fit$std.error))
})

test_that("a target must be valid", {
  p <- small_panel(seed = 43)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  expect_error(trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none",
                           target = last_period_off(p$W)), "only for `method = \"local\"`")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, se = "none", target = matrix(FALSE, nrow(p$W), ncol(p$W))),
               "no treated cell with an outcome")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, se = "none", target = matrix(TRUE, 2, 2)),
               "dimensions of `Y`")
})

test_that("trop() reads target and strata from columns of the data", {
  p <- small_panel(seed = 47, n_units = 16)
  d <- p$data
  d$in_target <- d$time < max(d$time)
  d$wave <- ifelse(d$unit %% 2 == 0, "a", "b")
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop(d, "y", "d", "unit", "time", lambda = l, target = "in_target", strata = "wave",
              n_boot = 10, seed = 1, control = tight_control())
  tgt <- matrix(d$in_target, nrow(p$W), ncol(p$W))
  ref <- trop_matrix(p$Y, p$W, lambda = l, target = tgt, strata = rep(c("b", "a"), 8), n_boot = 10,
                     seed = 1, control = tight_control())
  expect_equal(fit$estimate, ref$estimate, tolerance = 1e-10)
  expect_equal(fit$boot$draws, ref$boot$draws, tolerance = 1e-10)
  expect_true(fit$boot$strata)

  d$wave[1L] <- "c"
  expect_error(trop(d, "y", "d", "unit", "time", lambda = l, strata = "wave", n_boot = 5),
               "one value for each unit")
  expect_error(trop(p$data, "y", "d", "unit", "time", lambda = l, target = "nope"),
               "name of one column")
  d2 <- p$data
  d2$bad <- 2
  expect_error(trop(d2, "y", "d", "unit", "time", lambda = l, target = "bad"),
               "only TRUE and FALSE")
})

test_that("trop_boot_rows gives the draws of a fit with the same seed", {
  p <- small_panel(seed = 51)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop_matrix(p$Y, p$W, lambda = l, n_boot = 12, seed = 3, control = small_control())
  expect_equal(fit$boot$rows, trop_boot_rows(p$W, 12, seed = 3))
  set.seed(7)
  before <- .Random.seed
  trop_boot_rows(p$W, 5, seed = 1)
  expect_identical(.Random.seed, before)
})

test_that("strata keep the number of treated and untreated units in each stratum", {
  p <- small_panel(seed = 53, n_units = 16, n_treated = 6)
  N <- nrow(p$W)
  strata <- rep(c("a", "b", "c"), length.out = N)
  treated <- which(rowSums(p$W) > 0)
  rows <- trop_boot_rows(p$W, 25, strata = strata, seed = 2)
  expect_equal(dim(rows), c(25L, N))
  for (b in seq_len(nrow(rows))) {
    r <- rows[b, ]
    expect_equal(sum(r %in% treated), length(treated))
    expect_equal(c(table(strata[r], r %in% treated)),
                 c(table(strata, seq_len(N) %in% treated)))
  }
  rc <- trop_boot_rows(p$W, 10, strata = strata, resample = "controls", seed = 2)
  expect_true(all(apply(rc[, seq_along(treated), drop = FALSE], 1, identical, treated)))
})

test_that("strata and boot_rows must be valid and cannot be used together", {
  p <- small_panel(seed = 55)
  N <- nrow(p$W)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  rows <- trop_boot_rows(p$W, 6, seed = 1)
  expect_error(trop_boot_rows(p$W, 5, strata = 1:3), "one value for each unit")
  expect_error(trop_boot_rows(p$W, 5, strata = c(NA, rep(1, N - 1))), "no `NA` values")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows, strata = rep(1, N)),
               "not both")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows[, -1]), "boot_rows")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows[1, , drop = FALSE]), "at least 2 rows")
  bad <- rows
  bad[1, 1] <- N + 1L
  expect_error(trop_matrix(p$Y, p$W, lambda = l, boot_rows = bad), "row numbers")
})

test_that("boot_rows gives draws that equal re-estimation on the rows", {
  p <- small_panel(seed = 57)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  rows <- trop_boot_rows(p$W, 5, seed = 4)
  fit <- trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows, control = tight_control())
  expect_identical(fit$boot$resample, "custom")
  expect_equal(fit$boot$rows, rows)
  for (b in 1:5) {
    r <- rows[b, ]
    direct <- suppressWarnings(trop_matrix(p$Y[r, ], p$W[r, ], lambda = l, se = "none",
                                           control = tight_control()))
    expect_equal(fit$boot$draws[b], direct$estimate, tolerance = 1e-8)
  }
  # Two fits with the same rows have a joint distribution.
  fit2 <- trop_matrix(p$Y * 2, p$W, lambda = l / c(1, 2, 1), boot_rows = rows,
                      control = tight_control())
  expect_equal(fit2$boot$rows, fit$boot$rows)
})

test_that("demeaning does not change the effects when the unit weights are 0", {
  p <- small_panel(seed = 59)
  l <- c(time = 0.3, unit = 0, nn = 1)
  a <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = tight_control())
  b <- trop_matrix(p$Y, p$W, lambda = l, se = "none", demean = TRUE, control = tight_control())
  expect_equal(b$tau, a$tau, tolerance = 1e-7)
  expect_true(b$demean)
})

test_that("with demeaning, a constant added to each unit does not change the fit", {
  p <- small_panel(seed = 61)
  l <- c(time = 0.3, unit = 0.6, nn = 1)
  shift <- seq_len(nrow(p$Y)) * 10
  a <- trop_matrix(p$Y, p$W, lambda = l, se = "none", demean = TRUE, control = tight_control())
  b <- trop_matrix(p$Y + shift, p$W, lambda = l, se = "none", demean = TRUE,
                   control = tight_control())
  expect_equal(b$estimate, a$estimate, tolerance = 1e-7)
  expect_error(trop_matrix(p$Y, p$W, lambda = l, se = "none", demean = NA), "TRUE or FALSE")
})

test_that("trop_grid_scaled follows the scale of the outcome", {
  p <- small_panel(seed = 63)
  g1 <- trop_grid_scaled(p$Y, p$W)
  g100 <- trop_grid_scaled(100 * p$Y, p$W)
  s1 <- attr(g1, "scales")
  s100 <- attr(g100, "scales")
  expect_equal(s100$unit / s1$unit, 100, tolerance = 1e-8)
  expect_equal(s100$nn / s1$nn, 100, tolerance = 1e-8)
  base <- trop_grid()
  expect_equal(g1$time, base$time)
  expect_equal(g1$unit, base$unit / s1$unit)
  fin <- is.finite(base$nn)
  expect_equal(g1$nn[fin], base$nn[fin] * s1$nn)
  expect_identical(g1$nn[!fin], base$nn[!fin])
  expect_equal(g100$unit, g1$unit / 100, tolerance = 1e-8)
  expect_error(trop_grid_scaled(p$Y, p$W, demean = "yes"), "TRUE or FALSE")
})

test_that("grid_boundary names the finite values at the edge of the grid", {
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 1), nn = c(0.5, Inf))
  edge <- fasttrop:::grid_boundary
  expect_equal(edge(c(time = 0.5, unit = 0, nn = Inf), g), "largest value of `time` (0.5)")
  expect_equal(edge(c(time = 0, unit = 1, nn = 0.5), g),
               c("largest value of `unit` (1)", "smallest value of `nn` (0.5)"))
  expect_length(edge(c(time = 0.2, unit = 0.5, nn = 2), g), 0L)
  g1 <- trop_grid(time = 0.2, unit = 0.5, nn = 0.5)
  expect_length(edge(c(time = 0.2, unit = 0.5, nn = 0.5), g1), 0L)
})

test_that("the cross-validation warns when it chooses the edge of the grid", {
  p <- small_panel(seed = 19, n_units = 10, n_periods = 8, n_post = 2)
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 1), nn = c(0.5, Inf))
  expect_warning(
    fit <- trop_matrix(p$Y, p$W, grid = g, cv_cells = 20, se = "none", seed = 2,
                       control = small_control()),
    "chose the largest value of `time`"
  )
  expect_true(length(fit$cv$boundary) > 0L)
  g1 <- trop_grid(time = 0.2, unit = 0.5, nn = 0.5)
  expect_no_warning(fit1 <- trop_matrix(p$Y, p$W, grid = g1, se = "none", control = small_control()))
  expect_length(fit1$cv$boundary, 0L)
})

test_that("grid = \"auto\" uses the scaled grid", {
  p <- small_panel(seed = 65)
  Y <- 100 * p$Y
  fit <- suppressWarnings(trop_matrix(Y, p$W, grid = "auto", cv_cells = 40, se = "none", seed = 1,
                                      control = small_control()))
  g <- trop_grid_scaled(Y, p$W)
  expect_true(all(fit$cv$trace$unit %in% g$unit))
  expect_true(all(fit$cv$trace$nn %in% g$nn))
  expect_error(trop_matrix(Y, p$W, grid = "wide", se = "none"), "auto")
})

test_that("the print methods show the target, the strata, and the version", {
  p <- small_panel(seed = 67)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  fit <- trop_matrix(p$Y, p$W, lambda = l, target = last_period_off(p$W), n_boot = 10, seed = 1,
                     strata = rep(c("a", "b"), 7), control = small_control())
  expect_output(print(fit), "Cells outside of `target`: ")
  expect_output(print(summary(fit)), "within strata")
  rows <- trop_boot_rows(p$W, 6, seed = 1)
  fit2 <- trop_matrix(p$Y, p$W, lambda = l, boot_rows = rows, control = small_control())
  expect_output(print(summary(fit2)), "custom")
  expect_false(is.null(fit$package))
})
