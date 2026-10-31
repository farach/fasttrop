# Simulations with staggered adoption and treatment effects that differ
# between cohorts and between periods. In a staggered design, the treated
# units start the treatment in different periods. Theorem 5.1 of the paper is
# for a block design. Thus, these simulations do not test the theory. They
# measure the bias, the error, and the coverage of intervals in the setting
# that Hull (2026) studies.
#
# Each panel has N = 60 units and T = 16 periods. Four cohorts of 12 units
# start the treatment in periods 5, 8, 11, and 14, and stay treated. The 12
# other units are never treated. Thus, every treated cell has untreated cells
# in its own row and in its own column. All estimators use the same 360
# treated cells. The script stops if one of them skips a cell.
#   Simulation C: the effect is 1 in every treated cell. The untreated
#                 outcomes have unit effects, time effects, and noise. A
#                 random draw gives the cohort of each unit. All estimators are
#                 unbiased, and the simulation is a control for D and E.
#   Simulation D: as C, but the effect is delta + 0.1 * e. The value e is the
#                 number of periods since the start. The value delta is 0.8,
#                 0.6, 0.4, and 0.2 for the four cohorts. The mean of the effect
#                 over the 360 treated cells is 1, and the effect is between
#                 0.2 and 1.9. The untreated outcomes have parallel trends.
#   Simulation E: as D, and the untreated outcomes also have the two factors
#                 of simulation A. The units with the largest loading on the
#                 trend start first. The units with the smallest loadings are
#                 never treated. Thus, the trends of the cohorts differ.
#
# Methods, on the same panels:
#   - twfe_regression: the coefficient of the treatment indicator in the
#     regression of the outcome on unit effects, period effects, and the
#     indicator. The regression uses all 960 cells. The standard error is
#     robust to correlation in the units (clusters), with the factor G/(G-1).
#   - twfe_imputation: lambda = (0, 0, Inf). The fit uses only untreated cells,
#     and the estimate is the mean of the 360 treated-cell effects.
#   - local: the local method, with tuning by the leave-one-out criterion on 80
#     sampled untreated cells.
# The global method is only for a block design. Thus, the script does not use
# it. The standard errors of the two TROP methods come from the stratified unit
# bootstrap of Algorithm 3. The tuning parameters stay fixed in the bootstrap.
# The truth is the mean of the effects in the treated cells. The script gives
# the coverage of the 95% normal intervals (estimate +/- 1.96 SE) and of the
# percentile intervals. The regression has no percentile interval.
#
# The script also writes the weights that the regression puts on the treated
# cells. A weight of 1 is the equal weight of the imputation.
#
# Usage, in the package folder:
#   Rscript inst/validation/staggered.R [reps] [threads] [design]
library(fasttrop)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
threads <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
design <- if (length(args) >= 3L) args[3L] else "E"
stopifnot(design %in% c("C", "D", "E"))

N <- 60L
TT <- 16L
start <- c(5L, 8L, 11L, 14L)
per <- 12L
delta <- c(0.8, 0.6, 0.4, 0.2)
slope <- 0.1
sigma <- 0.5
boot_local <- 100L
out_dir <- file.path("inst", "validation")
set.seed(20261031)
f1 <- seq(-1, 1, length.out = TT)
f2 <- stats::rnorm(TT)

# Treatment indicator W and periods since the start E, from the first treated
# period of each unit. The value NA is for a unit that is never treated.
treatment <- function(adopt) {
  W <- matrix(0L, N, TT)
  E <- matrix(0L, N, TT)
  for (i in which(!is.na(adopt))) {
    on <- adopt[i]:TT
    W[i, on] <- 1L
    E[i, on] <- on - adopt[i]
  }
  list(W = W, E = E)
}

effect <- function(W, E, adopt) {
  if (design == "C") return(W * 1)
  d <- delta[match(adopt, start)]
  d[is.na(d)] <- 0
  (outer(d, rep(1, TT)) + slope * E) * W
}

draw_panel <- function() {
  a <- stats::rnorm(N)
  b <- cumsum(stats::rnorm(TT, 0, 0.3))
  G <- matrix(stats::rnorm(2L * N), N, 2L)
  eps <- matrix(stats::rnorm(N * TT, 0, sigma), N, TT)
  Y0 <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + eps
  if (design == "E") Y0 <- Y0 + outer(G[, 1], f1) + outer(G[, 2], f2)
  score <- if (design == "E") G[, 1] + stats::rnorm(N) else stats::rnorm(N)
  adopt <- integer(N)
  adopt[order(score, decreasing = TRUE)] <- rep(c(start, NA), times = c(rep(per, 4L), N - 4L * per))
  tr <- treatment(adopt)
  tau <- effect(tr$W, tr$E, adopt)
  list(Y = Y0 + tau, W = tr$W, truth = mean(tau[tr$W == 1L]))
}

# Two-way demeaning of a balanced panel. The coefficient of the treatment
# indicator in the regression with unit effects and period effects is the
# result of the Frisch-Waugh-Lovell theorem.
demean2 <- function(X) X - rowMeans(X) - rep(colMeans(X), each = nrow(X)) + mean(X)

twfe_regression <- function(Y, W) {
  Dt <- demean2(W)
  est <- sum(Dt * Y) / sum(Dt * W)
  e <- demean2(Y) - est * Dt
  se <- sqrt(nrow(Y) / (nrow(Y) - 1) * sum(rowSums(Dt * e)^2)) / sum(Dt * W)
  c(estimate = est, se = se)
}

# Weights of the regression on the treated cells, for the cohorts of the
# simulations. The weights add to the number of treated cells.
weights_file <- function() {
  adopt <- rep(c(start, NA), times = c(rep(per, 4L), N - 4L * per))
  tr <- treatment(adopt)
  Dt <- demean2(tr$W)
  w <- Dt * tr$W / sum(Dt * tr$W) * sum(tr$W)
  rows <- lapply(seq_along(start), function(k) {
    i <- which(adopt == start[k])[1L]
    on <- start[k]:TT
    data.frame(cohort_start = start[k], periods_since_start = on - start[k], period = on,
               weight = w[i, on])
  })
  do.call(rbind, rows)
}
utils::write.csv(weights_file(), file.path(out_dir, "staggered-twfe-weights.csv"),
                 row.names = FALSE)

grid <- trop_grid(time = c(0, 0.25, 0.5, 1), unit = c(0, 0.5, 1, 2), nn = c(0.1, 0.5, 2, Inf))
ctl <- trop_control(threads = threads)
fit_trop <- function(p, s, which) {
  suppressWarnings(
    if (which == "twfe_imputation") {
      trop_matrix(p$Y, p$W, lambda = c(time = 0, unit = 0, nn = Inf), n_boot = 200L, seed = s,
                  control = ctl)
    } else {
      trop_matrix(p$Y, p$W, method = "local", grid = grid, cv_cells = 80L, n_boot = boot_local,
                  seed = s, control = ctl)
    }
  )
}

t_start <- Sys.time()
rows <- vector("list", reps * 3L)
n <- 0L
for (r in seq_len(reps)) {
  p <- draw_panel()
  n_cells <- sum(p$W)
  t0 <- Sys.time()
  reg <- twfe_regression(p$Y, p$W)
  n <- n + 1L
  rows[[n]] <- data.frame(rep = r, method = "twfe_regression", estimate = reg[["estimate"]],
                          se = reg[["se"]], lo = NA_real_, hi = NA_real_, truth = p$truth,
                          lambda_time = NA_real_, lambda_unit = NA_real_, lambda_nn = NA_real_,
                          seconds = as.double(difftime(Sys.time(), t0, units = "secs")),
                          cells = n_cells, n_nonconverged = 0L, n_boot_failed = 0L)
  for (m in c("twfe_imputation", "local")) {
    fit <- fit_trop(p, r, m)
    stopifnot(fit$n_skipped == 0L, fit$n_treated_cells == n_cells)
    pct <- unname(stats::quantile(fit$boot$draws, c(0.025, 0.975), type = 7, na.rm = TRUE))
    n <- n + 1L
    rows[[n]] <- data.frame(rep = r, method = m, estimate = fit$estimate, se = fit$std.error,
                            lo = pct[1], hi = pct[2], truth = p$truth,
                            lambda_time = fit$lambda[["time"]], lambda_unit = fit$lambda[["unit"]],
                            lambda_nn = fit$lambda[["nn"]], seconds = fit$timing$total_seconds,
                            cells = fit$n_treated_cells - fit$n_skipped,
                            n_nonconverged = fit$n_nonconverged, n_boot_failed = fit$boot$n_failed)
  }
  if (r %% 10L == 0L) {
    message(format(Sys.time(), "%H:%M:%S"), " replication ", r, " of ", reps)
  }
}
res <- do.call(rbind, rows)
res$design <- design
utils::write.csv(res, file.path(out_dir, sprintf("staggered-draws-%s.csv", design)),
                 row.names = FALSE)

summarize <- function(d) {
  z <- stats::qnorm(0.975)
  err <- d$estimate - d$truth
  cov_n <- mean(abs(err) <= z * d$se)
  cov_p <- mean(d$lo <= d$truth & d$truth <= d$hi)
  data.frame(design = d$design[1], method = d$method[1], reps = nrow(d), cells = d$cells[1],
             bias = mean(err), mc_se_bias = stats::sd(d$estimate) / sqrt(nrow(d)),
             sd = stats::sd(d$estimate), mean_se = mean(d$se), rmse = sqrt(mean(err^2)),
             coverage_normal = cov_n, mc_se_normal = sqrt(cov_n * (1 - cov_n) / nrow(d)),
             coverage_percentile = cov_p,
             mc_se_percentile = sqrt(cov_p * (1 - cov_p) / nrow(d)),
             median_seconds = stats::median(d$seconds))
}
summ <- do.call(rbind, lapply(split(res, res$method), summarize))
summ <- summ[match(c("twfe_regression", "twfe_imputation", "local"), summ$method), ]
row.names(summ) <- NULL
utils::write.csv(summ, file.path(out_dir, sprintf("staggered-summary-%s.csv", design)),
                 row.names = FALSE)
message("Total time: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)))
message("Cells not converged: ", sum(res$n_nonconverged), ". Bootstrap draws failed: ",
        sum(res$n_boot_failed), ".")
print(summ, digits = 3)
