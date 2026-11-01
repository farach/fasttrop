# Simulation with a known effect and a time-varying covariate (Section 6.2 of
# the paper).
#
# The untreated outcomes have unit effects, time effects, two factors, noise,
# and the term X phi with phi = 1. The covariate is X_it = z_i s_t + e_it. The
# shock s_t is 0 in the pre-periods and 1.5 in the last k periods. The exposure
# z_i is larger for the treated units. Thus, the pre-periods do not show the
# shock. The unit weights and the low-rank term cannot remove the bias that it
# makes. Only the observed covariate can remove this bias. The treatment
# effect is 1 in the last k periods.
#   N = 40, T = 16, 10 treated units, k = 4.
#
# Methods, each without covariates and with covariates:
#   - TWFE imputation (lambda = (0, 0, Inf))
#   - The local method, with tuning by the leave-one-out criterion on 80
#     sampled untreated cells
#   - The global method, with tuning by placebo cross-validation on 50 sets
# The suffix "_joint" or "_resid" gives the covariate method. For TWFE, the
# two covariate methods give the same estimate. Thus, the script uses only
# "twfe_joint". The standard errors come from the stratified unit bootstrap of
# Algorithm 3.
#
# Usage, in the package folder:
#   Rscript inst/validation/covariates.R [reps] [threads]
library(fasttrop)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
threads <- if (length(args) >= 2L) as.integer(args[2L]) else 8L

N <- 40L
TT <- 16L
N1 <- 10L
k <- 4L
boot_local <- 100L
tau <- 1
phi <- 1
sigma <- 0.5
set.seed(20261101)
f1 <- seq(-1, 1, length.out = TT)
f2 <- stats::rnorm(TT)
shock <- c(rep(0, TT - k), rep(1.5, k))

draw_panel <- function() {
  a <- stats::rnorm(N)
  b <- cumsum(stats::rnorm(TT, 0, 0.3))
  G <- matrix(stats::rnorm(2L * N), N, 2L)
  treated <- sample.int(N, N1)
  z <- stats::rnorm(N)
  z[treated] <- z[treated] + 1
  X <- outer(z, shock) + matrix(stats::rnorm(N * TT, 0, 0.5), N, TT)
  Y0 <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + outer(G[, 1], f1) + outer(G[, 2], f2) +
    phi * X + matrix(stats::rnorm(N * TT, 0, sigma), N, TT)
  W <- matrix(0L, N, TT)
  W[treated, (TT - k + 1L):TT] <- 1L
  list(Y = Y0 + tau * W, W = W, X = array(X, c(N, TT, 1L), dimnames = list(NULL, NULL, "x")))
}

grid <- trop_grid(time = c(0, 0.25, 0.5, 1), unit = c(0, 0.5, 1, 2), nn = c(0.1, 0.5, 2, Inf))
ctl <- trop_control(threads = threads)
twfe <- function(p, s, X = NULL, cm = "joint") {
  trop_matrix(p$Y, p$W, lambda = c(time = 0, unit = 0, nn = Inf), n_boot = 200L, seed = s,
              control = ctl, X = X, covariate_method = cm)
}
local <- function(p, s, X = NULL, cm = "joint") {
  trop_matrix(p$Y, p$W, method = "local", grid = grid, cv_cells = 80L, n_boot = boot_local,
              seed = s, control = ctl, X = X, covariate_method = cm)
}
global <- function(p, s, X = NULL, cm = "joint") {
  trop_matrix(p$Y, p$W, method = "global", grid = grid, n_placebo = 50L, n_boot = 200L,
              seed = s, control = ctl, X = X, covariate_method = cm)
}
methods <- list(
  twfe = function(p, s) twfe(p, s),
  twfe_joint = function(p, s) twfe(p, s, p$X),
  local = function(p, s) local(p, s),
  local_joint = function(p, s) local(p, s, p$X),
  local_resid = function(p, s) local(p, s, p$X, "residualize"),
  global = function(p, s) global(p, s),
  global_joint = function(p, s) global(p, s, p$X),
  global_resid = function(p, s) global(p, s, p$X, "residualize")
)

t_start <- Sys.time()
rows <- vector("list", reps * length(methods))
n <- 0L
for (r in seq_len(reps)) {
  p <- draw_panel()
  for (m in names(methods)) {
    fit <- suppressWarnings(methods[[m]](p, r))
    pct <- unname(stats::quantile(fit$boot$draws, c(0.025, 0.975), type = 7, na.rm = TRUE))
    phi_hat <- if (is.null(fit$covariates)) NA_real_ else unname(fit$covariates$coef[1L])
    n <- n + 1L
    rows[[n]] <- data.frame(rep = r, method = m, estimate = fit$estimate, se = fit$std.error,
                            lo = pct[1], hi = pct[2], phi = phi_hat,
                            lambda_time = fit$lambda[["time"]],
                            lambda_unit = fit$lambda[["unit"]], lambda_nn = fit$lambda[["nn"]],
                            seconds = fit$timing$total_seconds)
  }
  if (r %% 10L == 0L) {
    message(format(Sys.time(), "%H:%M:%S"), " replication ", r, " of ", reps)
  }
}
res <- do.call(rbind, rows)
out_dir <- file.path("inst", "validation")
utils::write.csv(res, file.path(out_dir, "covariates-draws.csv"), row.names = FALSE)

summarize <- function(d) {
  z <- stats::qnorm(0.975)
  cov_n <- mean(abs(d$estimate - tau) <= z * d$se)
  cov_p <- mean(d$lo <= tau & tau <= d$hi)
  data.frame(method = d$method[1], reps = nrow(d),
             bias = mean(d$estimate) - tau, mc_se_bias = stats::sd(d$estimate) / sqrt(nrow(d)),
             sd = stats::sd(d$estimate), mean_se = mean(d$se),
             rmse = sqrt(mean((d$estimate - tau)^2)),
             coverage_normal = cov_n, mc_se_normal = sqrt(cov_n * (1 - cov_n) / nrow(d)),
             coverage_percentile = cov_p, mc_se_percentile = sqrt(cov_p * (1 - cov_p) / nrow(d)),
             mean_phi = mean(d$phi), median_seconds = stats::median(d$seconds))
}
summ <- do.call(rbind, lapply(split(res, res$method), summarize))
summ <- summ[match(names(methods), summ$method), ]
row.names(summ) <- NULL
utils::write.csv(summ, file.path(out_dir, "covariates-summary.csv"), row.names = FALSE)
message("Total time: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)))
print(summ, digits = 3)
