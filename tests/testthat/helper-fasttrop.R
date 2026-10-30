# Helper functions for the tests of fasttrop.

read_fixture <- function(name) {
  df <- utils::read.csv(test_path("fixtures", paste0(name, ".csv")))
  pm <- fasttrop:::panel_matrices(df, "y", "d", "unit", "time")
  list(data = df, Y = pm$Y, W = pm$W)
}

fixture_lambdas <- function() {
  utils::read.csv(test_path("fixtures", "lambdas.csv"))
}

lambda_of <- function(row) c(time = row$time, unit = row$unit, nn = row$nn)

# Strict test: the largest absolute difference must not be more than tol.
expect_max_abs_diff <- function(actual, expected, tol, label) {
  expect_equal(length(actual), length(expected), label = paste(label, "(length)"))
  expect_lte(max(abs(actual - expected)), tol, label = paste(label, "max |difference|"))
}

# Strict solver options for comparisons with the reference implementations.
tight_control <- function(...) {
  trop_control(tol = 1e-10, max_iter = 100000L, threads = 2L, ...)
}

small_control <- function(...) trop_control(threads = 2L, ...)

# A small panel with a block design and a factor structure. It has 5 treated
# units. Thus, the bootstrap does not give the warning for few treated units.
small_panel <- function(seed = 1, n_units = 14, n_periods = 10, n_treated = 5, n_post = 3,
                        design = "block") {
  sim <- trop_simulate(n_units = n_units, n_periods = n_periods, n_treated = n_treated,
                       n_post = n_post, effect = 1, design = design, seed = seed)
  list(data = sim, Y = matrix(sim$y, n_units, n_periods), W = matrix(sim$d, n_units, n_periods))
}

# Evaluate expr and return the messages of all its warnings.
collect_warnings <- function(expr) {
  msgs <- character()
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  msgs
}

# Independent R calculation of the local weights (Eq. 3 of the paper). The time
# weight is exp(-lt |t - s|). The unit weight is exp(-lu * RMS gap). The RMS
# gap uses the periods in which both units are observed controls, but not the
# target period t.
local_weights_r <- function(Y, W, i, t, lt, lu) {
  ctrl <- (W == 0) & is.finite(Y)
  theta <- exp(-lt * abs(t - seq_len(ncol(Y))))
  omega <- vapply(seq_len(nrow(Y)), function(j) {
    if (j == i || lu == 0) return(1)
    s <- setdiff(which(ctrl[i, ] & ctrl[j, ]), t)
    if (length(s) == 0L) return(0)
    exp(-lu * sqrt(mean((Y[i, s] - Y[j, s])^2)))
  }, double(1))
  outer(omega, theta) * ctrl
}

# Weighted two-way fixed effects prediction for cell (i, t) with lm().
wtwfe_predict <- function(Y, w, i, t) {
  d <- data.frame(y = as.vector(Y), unit = factor(rep(seq_len(nrow(Y)), ncol(Y))),
                  time = factor(rep(seq_len(ncol(Y)), each = nrow(Y))), w = as.vector(w))
  sub <- droplevels(d[d$w > 0, ])
  fit <- stats::lm(y ~ unit + time, data = sub, weights = w)
  nd <- d[i + nrow(Y) * (t - 1L), ]
  nd$unit <- factor(as.character(nd$unit), levels = levels(sub$unit))
  nd$time <- factor(as.character(nd$time), levels = levels(sub$time))
  unname(stats::predict(fit, newdata = nd))
}
