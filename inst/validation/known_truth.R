# Condition C3 of the council: simulations with a known effect. The
# simulations use a block design and a factor model. This is the setting of
# Theorem 5.1 of the paper.
#
# The untreated outcomes have unit effects, time effects, two factors, and
# noise. The first factor is a fixed linear trend. The script draws the second
# factor one time and then keeps it fixed. The treatment effect is 1 in the
# last k periods.
#   Simulation A:  N = 40, T = 16, 10 treated units, k = 4. The treated units
#                  are a random sample. Thus, all estimators are unbiased.
#   Simulation A2: as A, but N = 80, T = 12, 30 treated units, and k = 3. This
#                  simulation examines the bootstrap coverage with more treated
#                  units. The theory for Algorithm 3 of the paper is for a
#                  large number of treated units.
#   Simulation B:  as A, but the trend loading selects the treated units.
#                  Thus, the treated units have a different trend, and two-way
#                  fixed effects (TWFE) has a bias. This simulation measures
#                  how much of this bias TROP removes.
#
# Methods:
#   - TWFE (lambda = (0, 0, Inf)) with the same unit bootstrap
#   - The local method, with tuning by the leave-one-out criterion on 80
#     sampled untreated cells
#   - The global method, with tuning by placebo cross-validation on 50 sets
# The standard errors come from the stratified unit bootstrap of Algorithm 3.
# The tuning parameters stay fixed in the bootstrap. The script gives the
# coverage of the 95% normal intervals (estimate +/- 1.96 SE) and of the
# percentile intervals.
#
# Usage, in the package folder:
#   Rscript inst/validation/known_truth.R [reps] [threads] [design]
library(fasttrop)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
threads <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
design <- if (length(args) >= 3L) args[3L] else "B"
stopifnot(design %in% c("A", "A2", "B"))

if (design == "A2") {
  N <- 80L
  TT <- 12L
  N1 <- 30L
  k <- 3L
} else {
  N <- 40L
  TT <- 16L
  N1 <- 10L
  k <- 4L
}
boot_local <- 100L
tau <- 1
sigma <- 0.5
set.seed(20261030)
f1 <- seq(-1, 1, length.out = TT)
f2 <- stats::rnorm(TT)

draw_panel <- function() {
  a <- stats::rnorm(N)
  b <- cumsum(stats::rnorm(TT, 0, 0.3))
  G <- matrix(stats::rnorm(2L * N), N, 2L)
  Y0 <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + outer(G[, 1], f1) + outer(G[, 2], f2) +
    matrix(stats::rnorm(N * TT, 0, sigma), N, TT)
  score <- if (design == "B") G[, 1] + stats::rnorm(N) else stats::rnorm(N)
  treated <- order(score, decreasing = TRUE)[seq_len(N1)]
  W <- matrix(0L, N, TT)
  W[treated, (TT - k + 1L):TT] <- 1L
  list(Y = Y0 + tau * W, W = W)
}

grid <- trop_grid(time = c(0, 0.25, 0.5, 1), unit = c(0, 0.5, 1, 2), nn = c(0.1, 0.5, 2, Inf))
ctl <- trop_control(threads = threads)
methods <- list(
  twfe = function(p, s) {
    trop_matrix(p$Y, p$W, lambda = c(time = 0, unit = 0, nn = Inf), n_boot = 200L, seed = s,
                control = ctl)
  },
  local = function(p, s) {
    trop_matrix(p$Y, p$W, method = "local", grid = grid, cv_cells = 80L, n_boot = boot_local,
                seed = s, control = ctl)
  },
  global = function(p, s) {
    trop_matrix(p$Y, p$W, method = "global", grid = grid, n_placebo = 50L, n_boot = 200L,
                seed = s, control = ctl)
  }
)

t_start <- Sys.time()
rows <- vector("list", reps * length(methods))
n <- 0L
for (r in seq_len(reps)) {
  p <- draw_panel()
  for (m in names(methods)) {
    fit <- suppressWarnings(methods[[m]](p, r))
    pct <- unname(stats::quantile(fit$boot$draws, c(0.025, 0.975), type = 7, na.rm = TRUE))
    n <- n + 1L
    rows[[n]] <- data.frame(rep = r, method = m, estimate = fit$estimate, se = fit$std.error,
                            lo = pct[1], hi = pct[2], lambda_time = fit$lambda[["time"]],
                            lambda_unit = fit$lambda[["unit"]], lambda_nn = fit$lambda[["nn"]],
                            seconds = fit$timing$total_seconds)
  }
  if (r %% 10L == 0L) {
    message(format(Sys.time(), "%H:%M:%S"), " replication ", r, " of ", reps)
  }
}
res <- do.call(rbind, rows)
res$design <- design
out_dir <- file.path("inst", "validation")
utils::write.csv(res, file.path(out_dir, sprintf("known-truth-draws-%s.csv", design)),
                 row.names = FALSE)

summarize <- function(d) {
  z <- stats::qnorm(0.975)
  cov_n <- mean(abs(d$estimate - tau) <= z * d$se)
  cov_p <- mean(d$lo <= tau & tau <= d$hi)
  data.frame(design = d$design[1], method = d$method[1], reps = nrow(d),
             bias = mean(d$estimate) - tau, mc_se_bias = stats::sd(d$estimate) / sqrt(nrow(d)),
             sd = stats::sd(d$estimate), mean_se = mean(d$se),
             rmse = sqrt(mean((d$estimate - tau)^2)),
             coverage_normal = cov_n, mc_se_normal = sqrt(cov_n * (1 - cov_n) / nrow(d)),
             coverage_percentile = cov_p, mc_se_percentile = sqrt(cov_p * (1 - cov_p) / nrow(d)),
             median_seconds = stats::median(d$seconds))
}
summ <- do.call(rbind, lapply(split(res, res$method), summarize))
summ <- summ[match(names(methods), summ$method), ]
row.names(summ) <- NULL
utils::write.csv(summ, file.path(out_dir, sprintf("known-truth-summary-%s.csv", design)),
                 row.names = FALSE)
message("Total time: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)))
print(summ, digits = 3)
