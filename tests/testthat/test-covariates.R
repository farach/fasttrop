cov_panel <- function(seed = 3, N = 30, TT = 12) {
  set.seed(seed)
  a <- rnorm(N)
  b <- rnorm(TT)
  x1 <- matrix(rnorm(N * TT), N, TT) + a
  x2 <- matrix(runif(N * TT), N, TT) + outer(rep(1, N), b)
  W <- matrix(0L, N, TT)
  W[1:6, (TT - 3):TT] <- 1L
  Y <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + 2 * x1 - 1.5 * x2 + 0.7 * W +
    matrix(rnorm(N * TT, sd = 0.3), N, TT)
  Y[N - 5, 3] <- NA
  X <- array(c(x1, x2), c(N, TT, 2), dimnames = list(NULL, NULL, c("x1", "x2")))
  list(Y = Y, W = W, X = X, N = N, TT = TT)
}

lm_reference <- function(p) {
  df <- data.frame(y = as.vector(p$Y), u = factor(rep(seq_len(p$N), p$TT)),
                   t = factor(rep(seq_len(p$TT), each = p$N)), x1 = as.vector(p$X[, , 1]),
                   x2 = as.vector(p$X[, , 2]), w = as.vector(p$W))
  m <- stats::lm(y ~ u + t + x1 + x2, data = df[df$w == 0 & !is.na(df$y), ])
  r <- df$y - stats::predict(m, df)
  list(att = mean(r[df$w == 1]), phi = unname(stats::coef(m)[c("x1", "x2")]))
}

test_that("TWFE imputation with covariates matches lm() for all methods", {
  p <- cov_panel()
  ref <- lm_reference(p)
  l <- c(time = 0, unit = 0, nn = Inf)
  for (meth in c("local", "global")) for (cm in c("joint", "residualize")) {
    f <- trop_matrix(p$Y, p$W, method = meth, lambda = l, se = "none", control = tight_control(),
                     X = p$X, covariate_method = cm)
    expect_equal(f$estimate, ref$att, tolerance = 1e-9)
    expect_equal(unname(f$covariates$coef), ref$phi, tolerance = 1e-9)
  }
})

test_that("residualize equals a manual two-step fit", {
  p <- cov_panel()
  ref <- lm_reference(p)
  # Y* = Y - Xt phi, where Xt is X without its unit and period effects on the
  # untreated cells.
  mask <- p$W == 0L & !is.na(p$Y)
  Ys <- p$Y
  for (k in 1:2) {
    fe <- within_fe(p$X[, , k], mask)
    Xt <- p$X[, , k] - attr(fe, "a") - rep(attr(fe, "b"), each = p$N)
    Ys <- Ys - ref$phi[k] * Xt
  }
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  for (meth in c("local", "global")) {
    a <- trop_matrix(p$Y, p$W, method = meth, lambda = l, se = "none", control = tight_control(),
                     X = p$X, covariate_method = "residualize")
    b <- trop_matrix(Ys, p$W, method = meth, lambda = l, se = "none", control = tight_control())
    expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
  }
})

test_that("unit and period constants in a covariate do not change the results", {
  p <- cov_panel()
  X2 <- p$X
  X2[, , 1] <- X2[, , 1] + outer(rnorm(p$N, sd = 50), rep(1, p$TT)) +
    outer(rep(1, p$N), rnorm(p$TT, sd = 50))
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  for (meth in c("local", "global")) for (cm in c("joint", "residualize")) {
    a <- trop_matrix(p$Y, p$W, method = meth, lambda = l, se = "none", control = tight_control(),
                     X = p$X, covariate_method = cm)
    b <- trop_matrix(p$Y, p$W, method = meth, lambda = l, se = "none", control = tight_control(),
                     X = X2, covariate_method = cm)
    expect_equal(a$estimate, b$estimate, tolerance = 1e-7)
    expect_equal(a$covariates$coef, b$covariates$coef, tolerance = 1e-7)
  }
})

test_that("a covariate with a large unit part and a small time-varying part is used", {
  p <- cov_panel()
  set.seed(7)
  x <- outer(rnorm(p$N, sd = 1e3), rep(1, p$TT)) + matrix(rnorm(p$N * p$TT, sd = 1e-3), p$N, p$TT)
  Y <- p$Y + 2 * x
  l <- c(time = 0.1, unit = 0.1, nn = Inf)
  for (meth in c("local", "global")) for (cm in c("joint", "residualize")) {
    f <- trop_matrix(Y, p$W, method = meth, lambda = l, se = "none", control = tight_control(),
                     X = x, covariate_method = cm)
    expect_equal(f$covariates$used, "x1")
    expect_true(f$covariates$coef[[1]] != 0)
  }
})

test_that("covariates that the fixed effects absorb do not change the estimate", {
  p <- cov_panel()
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  base <- trop_matrix(p$Y, p$W, method = "local", lambda = l, se = "none", control = tight_control())
  for (Xbad in list(array(0, c(p$N, p$TT, 1)), outer(rnorm(p$N), rep(1, p$TT)),
                    outer(rep(1, p$N), rnorm(p$TT)))) {
    expect_warning(
      f <- trop_matrix(p$Y, p$W, method = "local", lambda = l, se = "none",
                       control = tight_control(), X = Xbad),
      "does not use these covariates"
    )
    expect_equal(f$estimate, base$estimate, tolerance = 1e-12)
    expect_length(f$covariates$coef, 0L)
  }
  X3 <- array(c(p$X, 2 * p$X[, , 1] - p$X[, , 2]), c(p$N, p$TT, 3),
              dimnames = list(NULL, NULL, c("x1", "x2", "x3")))
  expect_warning(f3 <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none",
                                   control = tight_control(), X = X3), "x3")
  f2 <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none",
                    control = tight_control(), X = p$X)
  expect_equal(f3$estimate, f2$estimate, tolerance = 1e-10)
})

test_that("the coefficients do not depend on the units of the covariates", {
  p <- cov_panel()
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  Xs <- p$X
  Xs[, , 1] <- 1000 * Xs[, , 1] + 5
  a <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none", control = tight_control(),
                   X = p$X)
  b <- trop_matrix(p$Y, p$W, method = "global", lambda = l, se = "none", control = tight_control(),
                   X = Xs)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
  expect_equal(unname(b$covariates$coef[1]) * 1000, unname(a$covariates$coef[1]), tolerance = 1e-7)
})

test_that("trop() with covariate columns equals trop_matrix()", {
  p <- cov_panel()
  df <- data.frame(unit = rep(seq_len(p$N), p$TT), time = rep(seq_len(p$TT), each = p$N),
                   y = as.vector(p$Y), d = as.vector(p$W), x1 = as.vector(p$X[, , 1]),
                   x2 = as.vector(p$X[, , 2]))
  df <- df[!is.na(df$y), ]
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  a <- trop(df, "y", "d", "unit", "time", lambda = l, se = "none", control = tight_control(),
            covariates = c("x1", "x2"))
  b <- trop_matrix(p$Y, p$W, lambda = l, se = "none", control = tight_control(), X = p$X)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-12)
  expect_equal(a$covariates$coef, b$covariates$coef, tolerance = 1e-12)
  expect_output(print(summary(a)), "Covariate coefficients")
  df$x1[5] <- NA
  expect_error(trop(df, "y", "d", "unit", "time", lambda = l, se = "none", covariates = "x1"),
               "NA")
  expect_error(trop(df, "y", "d", "unit", "time", lambda = l, se = "none", covariates = "zz"),
               "not in")
  df$f <- "a"
  expect_error(trop(df, "y", "d", "unit", "time", lambda = l, se = "none", covariates = "f"),
               "numeric")
})

test_that("covariate fits are the same for all thread counts and both solvers", {
  p <- cov_panel()
  l <- c(time = 0.3, unit = 0.4, nn = 0.05)
  for (meth in c("local", "global")) {
    a <- trop_matrix(p$Y, p$W, method = meth, lambda = l, n_boot = 20, seed = 1, X = p$X,
                     control = trop_control(threads = 1, tol = 1e-9, max_iter = 5000))
    b <- trop_matrix(p$Y, p$W, method = meth, lambda = l, n_boot = 20, seed = 1, X = p$X,
                     control = trop_control(threads = 3, tol = 1e-9, max_iter = 5000))
    expect_identical(a$estimate, b$estimate)
    expect_identical(a$boot$draws, b$boot$draws)
    s <- trop_matrix(p$Y, p$W, method = meth, lambda = l, se = "none", X = p$X,
                     control = trop_control(threads = 1, tol = 1e-10, max_iter = 20000, solver = "block"))
    expect_equal(s$estimate, a$estimate, tolerance = 1e-5)
  }
})

test_that("cross-validation and the bootstrap run with covariates", {
  p <- cov_panel()
  g <- trop_grid(time = c(0, 0.5), unit = c(0, 0.5), nn = c(0.1, Inf))
  for (meth in c("local", "global")) for (cm in c("joint", "residualize")) {
    f <- suppressWarnings(trop_matrix(p$Y, p$W, method = meth, grid = g, cv_cells = 40,
                                      n_placebo = 10, n_boot = 20, seed = 2, X = p$X,
                                      covariate_method = cm, control = trop_control(threads = 2)))
    expect_true(is.finite(f$estimate))
    expect_true(is.finite(f$std.error))
    expect_equal(f$boot$n_failed, 0L)
  }
})
