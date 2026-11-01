# The time cost of covariates. The script fits the same panels with 0, 2, and
# 5 covariates, with fixed tuning parameters and no bootstrap. It writes
# inst/bench/covariate-timings.csv.
#
# Usage, in the package folder:
#   Rscript inst/bench/covariate_timings.R [threads] [reps]
library(fasttrop)
args <- commandArgs(trailingOnly = TRUE)
threads <- if (length(args) >= 1L) as.integer(args[1L]) else 8L
reps <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
ctl <- trop_control(threads = threads)
lam <- c(time = 0.3, unit = 0.8, nn = 0.05)

sizes <- list(c(N = 100L, TT = 30L, N1 = 20L, k = 6L), c(N = 400L, TT = 50L, N1 = 80L, k = 10L))
rows <- list()
for (s in sizes) {
  set.seed(1)
  N <- s[["N"]]
  TT <- s[["TT"]]
  X <- array(stats::rnorm(N * TT * 5L), c(N, TT, 5L),
             dimnames = list(NULL, NULL, paste0("x", 1:5)))
  Y <- outer(stats::rnorm(N), rep(1, TT)) + outer(rep(1, N), stats::rnorm(TT)) +
    X[, , 1] - 0.5 * X[, , 2] + matrix(stats::rnorm(N * TT), N, TT)
  W <- matrix(0L, N, TT)
  W[seq_len(s[["N1"]]), (TT - s[["k"]] + 1L):TT] <- 1L
  Y <- Y + W
  for (p in c(0L, 2L, 5L)) {
    Xp <- if (p == 0L) NULL else X[, , seq_len(p), drop = FALSE]
    for (m in c("local", "global")) {
      reps_m <- if (m == "local" && N >= 400L) 1L else reps
      t <- vapply(seq_len(reps_m), function(r) {
        f <- trop_matrix(Y, W, method = m, lambda = lam, se = "none", control = ctl, X = Xp)
        c(f$timing$total_seconds, f$timing$kernel_seconds)
      }, double(2))
      rows[[length(rows) + 1L]] <- data.frame(
        n_units = N, n_periods = TT, n_treated_cells = sum(W), covariates = p, method = m,
        threads = threads, reps = reps_m, median_seconds = stats::median(t[1, ]),
        median_kernel_seconds = stats::median(t[2, ]))
      message(sprintf("N=%d T=%d p=%d %-6s %.3fs", N, TT, p, m, stats::median(t[1, ])))
    }
  }
}
res <- do.call(rbind, rows)
base <- res[res$covariates == 0L, c("n_units", "method", "median_seconds")]
names(base)[3] <- "base_seconds"
res <- merge(res, base, by = c("n_units", "method"), sort = FALSE)
res$ratio_to_no_covariates <- res$median_seconds / res$base_seconds
res$base_seconds <- NULL
res <- res[order(res$n_units, res$method, res$covariates), ]
utils::write.csv(res, file.path("inst", "bench", "covariate-timings.csv"), row.names = FALSE)
print(res, digits = 3)
