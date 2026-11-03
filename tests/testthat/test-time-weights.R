test_that("local time weights retain direct-exponential underflow boundaries", {
  for (distance in c(1490L, 1491L, 2000L)) {
    TT <- distance + 1L
    Y <- matrix(NA_real_, 3L, TT)
    Y[cbind(c(1, 2, 2, 3, 1), c(1, 1, TT, TT, TT))] <- c(1, 2, 3, 4, 5)
    W <- matrix(0L, 3L, TT)
    W[1, TT] <- 1L
    lambda <- c(time = 0.5, unit = 0, nn = Inf)
    positive <- exp(-lambda[["time"]] * distance) > 0
    for (reverse in c(FALSE, TRUE)) {
      columns <- if (reverse) rev(seq_len(TT)) else seq_len(TT)
      inp <- fasttrop:::as_panel_inputs(Y[, columns], W[, columns])
      target <- if (reverse) 1L else TT
      cell <- as.integer(3L * (target - 1L))
      for (threads in c(1L, 2L)) {
        ctl <- unclass(trop_control(threads = threads, max_iter = 10L))
        est <- .Call(fasttrop:::C_local_estimate, inp$Y, inp$W, inp$obs,
                     NULL, lambda, ctl)
        expect_identical(est$n_used, as.integer(positive))
        expect_identical(est$n_skipped, as.integer(!positive))
        if (!positive) expect_true(is.nan(est$tau[1, target]))

        Wcv <- matrix(0L, 3L, TT)
        cv <- .Call(fasttrop:::C_local_loocv, inp$Y, Wcv, inp$obs,
                    NULL, lambda, cell, ctl)
        expect_identical(cv$n_valid, as.integer(positive))
        if (!positive) expect_identical(cv$Q, Inf)

        search <- .Call(fasttrop:::C_local_search, inp$Y, Wcv, inp$obs,
                        NULL, c(0.5, 0), 0, Inf, cell, 2L, ctl)
        at_decay <- search$trace[, 1L] == 0.5
        expect_true(any(at_decay))
        expect_equal(search$trace[at_decay, 4L], cv$Q, tolerance = 0)

        rows <- matrix(rep(0:2, each = 2L), 2L, 3L)
        boot <- .Call(fasttrop:::C_bootstrap, inp$Y, inp$W, inp$obs,
                      NULL, FALSE, 0L, 0L, lambda, rows, ctl)
        expect_equal(boot$draws, rep(est$att, 2L), tolerance = 0)
      }
    }
  }
})

test_that("cached CV scores equal separate fits and do not depend on thread count", {
  p <- small_panel(seed = 147, n_units = 8, n_periods = 9, n_treated = 2)
  p$Y[3, 2] <- NA_real_
  inp <- fasttrop:::as_panel_inputs(p$Y, p$W)
  cells <- as.integer(which(inp$W == 0L & inp$obs == 1L)[c(1, 8, 20, 35, 42)] - 1L)
  previous <- NULL
  for (threads in c(1L, 2L)) {
    ctl <- unclass(trop_control(threads = threads, tol = 1e-9))
    # Zero is absent from the time grid but is needed in the first stage.
    search <- .Call(fasttrop:::C_local_search, inp$Y, inp$W, inp$obs, NULL,
                    c(0.2, 0.7, 0.2), c(0, 0.8), c(0.5, Inf), cells, 3L, ctl)
    for (k in seq_len(nrow(search$trace))) {
      row <- search$trace[k, ]
      cv <- .Call(fasttrop:::C_local_loocv, inp$Y, inp$W, inp$obs,
                  NULL, unname(row[1:3]), cells, ctl)
      expect_equal(unname(row[4L]), cv$Q, tolerance = 0)
    }
    if (!is.null(previous)) {
      expect_identical(search$lambda, previous$lambda)
      expect_identical(search$score, previous$score)
      expect_identical(search$trace, previous$trace)
    }
    previous <- search
  }
})

test_that("local estimates agree with independently calculated exponential weights", {
  p <- small_panel(seed = 97, n_units = 7, n_periods = 8, n_treated = 2)
  p$Y[4, 3] <- NA_real_
  inp <- fasttrop:::as_panel_inputs(p$Y, p$W)
  ctl <- trop_control(threads = 1L, warm = FALSE, tol = 1e-10)
  for (lt in c(0, 0.3, 1, 1000)) {
    lambda <- c(time = lt, unit = 0.8, nn = 0.5)
    est <- .Call(fasttrop:::C_local_estimate, inp$Y, inp$W, inp$obs,
                 NULL, lambda, unclass(ctl))
    for (cell in which(p$W == 1L)) {
      i <- (cell - 1L) %% nrow(p$Y) + 1L
      t <- (cell - 1L) %/% nrow(p$Y) + 1L
      w <- local_weights_r(p$Y, p$W, i, t, lt, lambda[["unit"]])
      if (lt == 1000) {
        expect_true(is.nan(est$tau[i, t]))
      } else {
        Y <- p$Y
        Y[!is.finite(Y)] <- 0
        direct <- fasttrop:::trop_fit_weighted(Y, w, lambda[["nn"]], ctl)
        tau <- Y[i, t] - direct$alpha[i] - direct$beta[t] - direct$L[i, t]
        expect_equal(est$tau[i, t], tau, tolerance = 1e-8)
      }
    }
  }
})
