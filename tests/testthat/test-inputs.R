test_that("matrix inputs are validated with clear errors", {
  p <- small_panel()
  l <- c(time = 0, unit = 0, nn = Inf)
  expect_error(trop_matrix(as.data.frame(p$Y), p$W), "numeric matrix")
  expect_error(trop_matrix(p$Y, p$W[, -1]), "same dimensions")
  W2 <- p$W
  W2[1, 1] <- 2L
  expect_error(trop_matrix(p$Y, W2), "only 0 and 1")
  expect_error(trop_matrix(p$Y[1:2, ], p$W[1:2, ]), "at least 3 units")
  expect_error(trop_matrix(p$Y, p$W * 0L, lambda = l, se = "none"), "No treated")
  expect_error(trop_matrix(p$Y, p$W, lambda = c(1, 2)), "lambda")
  expect_error(trop_matrix(p$Y, p$W, lambda = c(time = -1, unit = 0, nn = 1)), "non-negative")
  expect_error(trop_matrix(p$Y, p$W, lambda = c(time = Inf, unit = 0, nn = 1)), "finite")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, conf_level = 1.5), "conf_level")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, ci = "bca"), "should be one of")
  expect_error(trop_matrix(p$Y, p$W, lambda = l, n_boot = 0), "n_boot")
  expect_error(trop_matrix(p$Y, p$W, cv_cells = -3, se = "none"), "cv_cells")
})

test_that("controls and grids are validated", {
  expect_error(trop_control(tol = 0), "tol")
  expect_error(trop_control(threads = 0), "threads")
  expect_error(trop_control(max_iter = 1.5), "max_iter")
  expect_error(trop_control(warm = NA), "warm")
  expect_error(trop_control(throttle = NA), "throttle")
  expect_error(trop_control(solver = "newton"), "should be one of")
  expect_error(trop_grid(time = -1), "time")
  expect_error(trop_grid(unit = Inf), "unit")
  expect_error(trop_grid(nn = NA), "nn")
  expect_s3_class(trop_control(threads = 1), "trop_control")
  expect_equal(trop_grid(nn = c(1, Inf))$nn, c(1, Inf))
})

test_that("data frame inputs are validated", {
  d <- small_panel()$data
  expect_error(trop(d, "nope", "d", "unit", "time"), "not found")
  d2 <- d
  d2$y <- as.character(d2$y)
  expect_error(trop(d2, "y", "d", "unit", "time"), "numeric")
  expect_error(trop(rbind(d, d[1, ]), "y", "d", "unit", "time"), "duplicated")
  d4 <- d
  d4$d[1] <- NA
  expect_error(trop(d4, "y", "d", "unit", "time"), "0/1")
  expect_error(trop(as.list(d), "y", "d", "unit", "time"), "data frame")
})

test_that("unbalanced panels are completed with missing outcomes", {
  d <- small_panel(seed = 2)$data
  d <- d[-c(3, 40), ] # drop two untreated rows
  fit <- trop(d, "y", "d", "unit", "time", lambda = c(time = 0, unit = 0, nn = Inf),
              se = "none", control = small_control())
  expect_true(is.finite(fit$estimate))
  expect_equal(fit$n_treated_cells, sum(d$d))
})

test_that("design checks warn or stop as documented", {
  l <- c(time = 0, unit = 0, nn = Inf)
  st <- small_panel(seed = 8, design = "staggered", n_post = 4)
  expect_equal(fasttrop:::trop_design(st$W, st$W * 0L + 1L)$type, "staggered")
  expect_warning(trop_matrix(st$Y, st$W, lambda = l, se = "none", control = small_control()),
                 "block assignment only")
  expect_error(suppressWarnings(trop_matrix(st$Y, st$W, method = "global", lambda = l,
                                            se = "none")),
               "block design")

  gen <- small_panel(seed = 9)
  i <- which(rowSums(gen$W) > 0)[1]
  gen$W[i, ncol(gen$W)] <- 0L
  warns <- collect_warnings(trop_matrix(gen$Y, gen$W, lambda = l, se = "none",
                                        control = small_control()))
  expect_true(any(grepl("switches off", warns)))
  expect_true(any(grepl("block assignment only", warns)))

  few <- small_panel(seed = 10, n_treated = 2)
  expect_warning(trop_matrix(few$Y, few$W, lambda = l, n_boot = 5, control = small_control()),
                 "growing number of treated units")
  one_control <- small_panel(seed = 12, n_units = 6, n_treated = 5)
  expect_error(trop_matrix(one_control$Y, one_control$W, method = "global", lambda = l,
                           se = "none"),
               "at least 2 never-treated")
})
