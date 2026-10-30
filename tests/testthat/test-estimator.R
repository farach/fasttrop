test_that("lambda = (0, 0, Inf) reproduces two-way fixed effects imputation", {
  p <- small_panel(seed = 3)
  Y <- p$Y
  W <- p$W
  Y[2, 3] <- NA # a missing untreated outcome
  d <- data.frame(y = as.vector(Y), unit = factor(rep(seq_len(nrow(Y)), ncol(Y))),
                  time = factor(rep(seq_len(ncol(Y)), each = nrow(Y))), w = as.vector(W))
  ols <- stats::lm(y ~ unit + time, data = d[d$w == 0 & is.finite(d$y), ])
  tr <- which(d$w == 1)
  tau <- d$y[tr] - stats::predict(ols, newdata = d[tr, ])
  for (method in c("local", "global")) {
    fit <- trop_matrix(Y, W, method = method, lambda = c(time = 0, unit = 0, nn = Inf),
                       se = "none", control = tight_control())
    expect_equal(fit$estimate, mean(tau), tolerance = 1e-8, label = paste(method, "ATT"))
    expect_equal(fit$tau[W == 1], unname(tau), tolerance = 1e-8, label = paste(method, "tau"))
  }
})

test_that("local weights follow the paper's definition (independent R check)", {
  p <- small_panel(seed = 5)
  Y <- p$Y
  Y[1, 2] <- NA
  lt <- 0.4
  lu <- 0.7
  fit <- trop_matrix(Y, p$W, lambda = c(time = lt, unit = lu, nn = Inf), se = "none",
                     control = tight_control())
  cells <- which(p$W == 1, arr.ind = TRUE)
  for (k in seq_len(nrow(cells))) {
    i <- cells[k, 1]
    t <- cells[k, 2]
    w <- local_weights_r(Y, p$W, i, t, lt, lu)
    expect_equal(fit$tau[i, t], Y[i, t] - wtwfe_predict(Y, w, i, t), tolerance = 1e-8,
                 label = sprintf("tau[%d, %d]", i, t))
  }
})

test_that("local and global coincide when unit and time weights are uniform", {
  p <- small_panel(seed = 43)
  l <- c(time = 0, unit = 0, nn = 0.5)
  a <- trop_matrix(p$Y, p$W, method = "local", lambda = l, se = "none", control = tight_control())
  b <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none", control = tight_control())
  expect_equal(a$estimate, b$estimate, tolerance = 1e-7)
})

test_that("the prox and block solvers agree at convergence", {
  p <- small_panel(seed = 7, n_units = 10, n_periods = 8, n_post = 2)
  l <- c(time = 0.3, unit = 0.5, nn = 1)
  a <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = tight_control(solver = "prox"))
  b <- trop_matrix(p$Y, p$W, lambda = l, se = "none",
                   control = trop_control(tol = 1e-10, max_iter = 50000L, threads = 2L,
                                          solver = "block"))
  expect_equal(b$n_nonconverged, 0L)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-6)
  expect_equal(a$tau, b$tau, tolerance = 1e-6)
})

test_that("results do not depend on the thread count or the throttle setting", {
  p <- small_panel(seed = 11)
  l <- c(time = 0.3, unit = 0.5, nn = 0.5)
  fit <- function(...) trop_matrix(p$Y, p$W, lambda = l, n_boot = 12, seed = 3,
                                   control = trop_control(...))
  f1 <- fit(threads = 1L)
  f2 <- fit(threads = 2L)
  f3 <- fit(threads = 2L, throttle = TRUE)
  expect_identical(f1$estimate, f2$estimate)
  expect_identical(f1$boot$draws, f2$boot$draws)
  expect_identical(f2$boot$draws, f3$boot$draws)
})

test_that("cells whose effects are not identified are skipped with a warning", {
  Y <- matrix(sin(1:20), 5, 4)
  W <- matrix(0L, 5, 4)
  W[1, 2:4] <- 1L # unit 1 has a single untreated period
  W[2:5, 4] <- 1L # period 4 has no untreated cell
  warns <- collect_warnings(
    fit <- trop_matrix(Y, W, lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                       control = small_control())
  )
  expect_true(any(grepl("skipped", warns)))
  expect_gt(fit$n_skipped, 0L)
  expect_true(all(is.na(fit$tau[, 4])))
})
