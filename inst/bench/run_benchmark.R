# Parity benchmark: fasttrop's R entry points versus a standalone C driver.
#
# Both paths execute the same compiled core (src/trop_core.c), built with the
# same compiler, flags, BLAS/LAPACK, and thread count, on the same inputs. The
# plan (panel, sampled CV cells, placebo sets, bootstrap rows) is drawn once in
# R and serialized for the driver, so timing differences come only from the R
# layer. Paths are interleaved in rounds so drift in machine load affects all
# of them alike.
#
# Usage (from the package root, with fasttrop installed):
#   Rscript inst/bench/run_benchmark.R [rounds] [threads] [inner]
#     rounds  number of interleaved rounds (default 3)
#     threads comma-separated thread counts (default "1,4")
#     inner   repetitions per path and round; the round keeps the minimum (default 2)
args <- commandArgs(trailingOnly = TRUE)
rounds <- if (length(args) >= 1L) as.integer(args[1L]) else 3L
thread_set <- if (length(args) >= 2L) as.integer(strsplit(args[2L], ",")[[1L]]) else c(1L, 4L)
inner <- if (length(args) >= 3L) as.integer(args[3L]) else 2L
library(fasttrop)
ns <- asNamespace("fasttrop")
log_msg <- function(...) {
  message(format(Sys.time(), "%H:%M:%S "), ...)
  flush(stderr())
}
bench_dir <- file.path(tempdir(), "fasttrop-bench")
dir.create(bench_dir, showWarnings = FALSE, recursive = TRUE)

# Build the driver with R's own toolchain configuration ------------------
rcmd <- function(x) {
  r <- file.path(R.home("bin"), "R")
  trimws(system2(r, c("CMD", "config", x), stdout = TRUE))
}
src_dir <- normalizePath(file.path(getwd(), "src"))
driver_c <- normalizePath(file.path(getwd(), "inst", "bench", "trop_driver.c"))
os_c <- normalizePath(file.path(src_dir, "trop_os.c"))
exe <- file.path(bench_dir, if (.Platform$OS.type == "windows") "trop_driver.exe" else "trop_driver")
make_conf <- file.path(R.home("etc"), Sys.getenv("R_ARCH"), "Makeconf")
omp <- sub("^SHLIB_OPENMP_CFLAGS\\s*=\\s*", "",
           grep("^SHLIB_OPENMP_CFLAGS\\s*=", readLines(make_conf), value = TRUE)[1L])
cc <- rcmd("CC")
cflags <- rcmd("CFLAGS")
libs <- paste(rcmd("LAPACK_LIBS"), rcmd("BLAS_LIBS"), rcmd("FLIBS"))
# R CMD INSTALL adds -DNDEBUG; do the same so both builds are identical.
cmd <- paste(cc, "-DNDEBUG", cflags, omp, sprintf('-I"%s"', R.home("include")),
             sprintf('-I"%s"', src_dir), sprintf('"%s"', driver_c), sprintf('"%s"', os_c),
             "-o", sprintf('"%s"', exe), libs, omp, "-lm")
log_msg("Compiling driver: ", cmd)
if (system(cmd) != 0L) stop("Driver compilation failed.")
if (.Platform$OS.type == "windows") {
  Sys.setenv(PATH = paste(R.home("bin"), Sys.getenv("PATH"), sep = .Platform$path.sep))
}

# Serialize a plan for the driver -----------------------------------------
write_plan <- function(plan, path) {
  con <- file(path, "wb")
  on.exit(close(con))
  wi <- function(x) writeBin(as.integer(x), con, size = 4L, endian = "little")
  wd <- function(x) writeBin(as.double(x), con, size = 8L, endian = "little")
  p <- plan$p
  ctrl <- plan$ctrl
  wi(0x54524F50L)
  wi(c(nrow(p$Y), ncol(p$Y), if (plan$method == "local") 0L else 1L,
       if (plan$method == "global") plan$design$k else 0L))
  wi(c(ctrl$threads, ctrl$max_iter, ctrl$max_inner, ctrl$solver, ctrl$warm, plan$max_cycles,
       ctrl$throttle))
  wd(ctrl$tol)
  wi(c(length(plan$grid$time), length(plan$grid$unit), length(plan$grid$nn)))
  wd(plan$grid$time); wd(plan$grid$unit); wd(plan$grid$nn)
  wi(as.integer(!is.null(plan$lambda)))
  wd(if (is.null(plan$lambda)) c(0, 0, 0) else plan$lambda)
  wi(length(plan$cells)); if (length(plan$cells)) wi(plan$cells)
  if (is.null(plan$Yc)) {
    wi(0L)
  } else {
    wi(nrow(plan$Yc)); wd(as.vector(plan$Yc)); wi(as.vector(plan$oc))
  }
  if (is.null(plan$sets)) wi(c(0L, 0L)) else { wi(dim(plan$sets)); wi(as.vector(t(plan$sets))) }
  wd(if (is.null(plan$init)) c(0, 0, 0) else plan$init)
  wi(length(plan$design$treated)); wi(plan$design$treated - 1L)
  if (is.null(plan$rows)) wi(0L) else { wi(nrow(plan$rows)); wi(as.vector(t(plan$rows))) }
  wd(as.vector(p$Y)); wi(as.vector(p$W)); wi(as.vector(p$obs))
}

run_driver <- function(plan_path, out_path, reps, B) {
  status <- system2(exe, c(shQuote(plan_path), reps, shQuote(out_path)), stdout = FALSE)
  if (!identical(as.integer(status), 0L)) stop("Driver failed with status ", status, ".")
  con <- file(out_path, "rb")
  on.exit(close(con))
  head <- readBin(con, "double", 8L)
  draws <- readBin(con, "double", B)
  list(lambda = head[1:3], att = head[4], total = head[5], search = head[6],
       estimate = head[7], boot = head[8], draws = draws)
}

# proc.time() has a 10 ms resolution on Windows; use a high-resolution clock.
if (!requireNamespace("bench", quietly = TRUE)) stop("The benchmark needs the 'bench' package.")
min_time <- function(f, n) {
  min(vapply(seq_len(n), function(i) {
    t0 <- bench::hires_time()
    f()
    as.double(bench::hires_time() - t0)
  }, double(1)))
}

# Benchmark configurations ------------------------------------------------
grid_small <- trop_grid(time = c(0, 0.5, 1), unit = c(0, 0.5, 1), nn = c(0.1, 1, Inf))
# `inner` is raised for the millisecond-scale configuration so the minimum is stable.
configs <- list(
  list(id = "A", name = "local, fixed lambda, no SE, 20x12", method = "local", N = 20, TT = 12,
       ntr = 4, post = 3, lambda = c(time = 0.3, unit = 0.8, nn = 0.5), cv_cells = NULL,
       se = "none", n_boot = 0L, inner = 25L),
  list(id = "B", name = "local, fixed lambda + bootstrap (B=100), 30x20", method = "local",
       N = 30, TT = 20, ntr = 6, post = 4, lambda = c(time = 0.3, unit = 0.8, nn = 0.5),
       cv_cells = NULL, se = "bootstrap", n_boot = 100L),
  list(id = "C", name = "local, LOOCV (60 cells) + bootstrap (B=50), 30x20", method = "local",
       N = 30, TT = 20, ntr = 6, post = 4, lambda = NULL, cv_cells = 60L, se = "bootstrap",
       n_boot = 50L),
  list(id = "D", name = "global, placebo CV (50 sets) + bootstrap (B=200), 50x30",
       method = "global", N = 50, TT = 30, ntr = 8, post = 6, lambda = NULL, cv_cells = NULL,
       se = "bootstrap", n_boot = 200L)
)

out_csv <- file.path(getwd(), "inst", "bench", "results.csv")
if (file.exists(out_csv)) file.remove(out_csv)
set.seed(20261030) # only for the order of paths within rounds
reference <- list()
results <- list()
paired <- list()
for (cf in configs) {
  sim <- trop_simulate(n_units = cf$N, n_periods = cf$TT, n_treated = cf$ntr, n_post = cf$post,
                       effect = 1, seed = 2026)
  Y <- matrix(sim$y, cf$N, cf$TT)
  W <- matrix(sim$d, cf$N, cf$TT)
  for (th in thread_set) {
    ctl <- trop_control(threads = th)
    args_fit <- list(Y = Y, W = W, method = cf$method, lambda = cf$lambda, grid = grid_small,
                     cv_cells = cf$cv_cells, n_placebo = 50L, se = cf$se,
                     n_boot = max(cf$n_boot, 1L), control = ctl, seed = 7)
    plan <- do.call(ns$trop_prepare, args_fit)
    plan_path <- file.path(bench_dir, "plan.bin")
    out_path <- file.path(bench_dir, "out.bin")
    write_plan(plan, plan_path)
    B <- if (is.null(plan$rows)) 0L else nrow(plan$rows)

    # Correctness first: the driver and R must agree exactly.
    drv <- run_driver(plan_path, out_path, 1L, B)
    fit <- suppressWarnings(ns$trop_run(plan))
    same_lambda <- isTRUE(all.equal(unname(fit$lambda), drv$lambda, tolerance = 0))
    att_diff <- abs(fit$estimate - drv$att)
    draw_diff <- if (B > 0L) max(abs(fit$boot$draws - drv$draws)) else 0
    key <- cf$id
    if (is.null(reference[[key]])) {
      reference[[key]] <- fit
      thread_att_diff <- 0
      thread_draw_diff <- 0
    } else {
      thread_att_diff <- abs(fit$estimate - reference[[key]]$estimate)
      thread_draw_diff <- if (B > 0L) max(abs(fit$boot$draws - reference[[key]]$boot$draws)) else 0
    }

    n_in <- if (is.null(cf$inner)) inner else cf$inner
    t_drv <- t_run <- t_full <- double(rounds)
    for (r in seq_len(rounds)) {
      # Random order within each round, so no path always runs first.
      for (path in sample(c("c", "run", "full"))) {
        if (path == "c") {
          t_drv[r] <- run_driver(plan_path, out_path, n_in, B)$total
        } else if (path == "run") {
          t_run[r] <- min_time(function() suppressWarnings(ns$trop_run(plan)), n_in)
        } else {
          t_full[r] <- min_time(function() suppressWarnings(do.call(trop_matrix, args_fit)), n_in)
        }
      }
      log_msg(sprintf("%s threads=%d round %d: C %.4fs  trop_run %.4fs  trop_matrix %.4fs",
                      cf$id, th, r, t_drv[r], t_run[r], t_full[r]))
    }
    row <- data.frame(
      config = cf$name, threads = th, rounds = rounds, inner = n_in,
      c_driver_s = stats::median(t_drv), r_trop_run_s = stats::median(t_run),
      r_trop_matrix_s = stats::median(t_full),
      overhead_run_ms = 1000 * (stats::median(t_run) - stats::median(t_drv)),
      overhead_run_pct = 100 * (stats::median(t_run) / stats::median(t_drv) - 1),
      overhead_full_ms = 1000 * (stats::median(t_full) - stats::median(t_drv)),
      overhead_full_pct = 100 * (stats::median(t_full) / stats::median(t_drv) - 1),
      c_driver_range_s = sprintf("%.4f-%.4f", min(t_drv), max(t_drv)),
      r_trop_run_range_s = sprintf("%.4f-%.4f", min(t_run), max(t_run)),
      identical_lambda = same_lambda, att_diff_vs_c = att_diff, draw_diff_vs_c = draw_diff,
      att_diff_vs_first_threads = thread_att_diff, draw_diff_vs_first_threads = thread_draw_diff
    )
    results[[length(results) + 1L]] <- row
    paired[[length(paired) + 1L]] <- data.frame(config = cf$id, threads = th, round = seq_len(rounds),
                                                c_s = t_drv, run_s = t_run, full_s = t_full)
    utils::write.csv(do.call(rbind, results), out_csv, row.names = FALSE)
    log_msg(sprintf("%s threads=%d done: overhead trop_run %+.2f%% (%+.1f ms), full %+.2f%%; lambda identical=%s, |att diff|=%.3g, max |draw diff|=%.3g, thread diff=%.3g",
                    cf$id, th, row$overhead_run_pct, row$overhead_run_ms, row$overhead_full_pct,
                    same_lambda, att_diff, draw_diff, max(thread_att_diff, thread_draw_diff)))
  }
}
log_msg("Saved ", out_csv)
print(do.call(rbind, results)[, c("config", "threads", "c_driver_s", "r_trop_run_s",
                                  "r_trop_matrix_s", "overhead_run_pct", "overhead_full_pct")])

# Pooled paired comparison for runs of at least one second: geometric mean of
# the per-round ratio R / C with a t-based 95% interval on the log scale.
pr <- do.call(rbind, paired)
utils::write.csv(pr, file.path(getwd(), "inst", "bench", "results-rounds.csv"), row.names = FALSE)
long <- pr[pr$c_s >= 1, ]
pooled <- do.call(rbind, lapply(c(trop_run = "run_s", trop_matrix = "full_s"), function(col) {
  lr <- log(long[[col]] / long$c_s)
  h <- stats::qt(0.975, length(lr) - 1L) * stats::sd(lr) / sqrt(length(lr))
  data.frame(n_pairs = length(lr), geo_mean_ratio = exp(mean(lr)),
             ci_low = exp(mean(lr) - h), ci_high = exp(mean(lr) + h))
}))
pooled$path <- rownames(pooled)
utils::write.csv(pooled, file.path(getwd(), "inst", "bench", "results-pooled.csv"), row.names = FALSE)
print(pooled, digits = 4)
