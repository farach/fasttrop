# Compare two installed package libraries on the same local CV searches.
# Usage: Rscript inst/bench/time_cache_timings.R BASE_LIB NEW_LIB [rounds]
# Each worker loads one library in a fresh R process. The C timer excludes
# process startup, package loading, panel preparation, and the warm-up run.
args <- commandArgs(trailingOnly = TRUE)

worker <- function(lib, output) {
  .libPaths(c(lib, .libPaths()))
  library(fasttrop)
  ns <- asNamespace("fasttrop")
  configs <- list(
    list(id = "30x20_cv60", N = 30L, TT = 20L, tr = 6L, post = 4L,
         cells = 60L, nn = c(0.1, 1, Inf), inner = 1L),
    list(id = "12x120_cv60_twfe", N = 12L, TT = 120L, tr = 3L, post = 12L,
         cells = 60L, nn = Inf, inner = 25L)
  )
  results <- list()
  for (cf in configs) {
    sim <- trop_simulate(n_units = cf$N, n_periods = cf$TT, n_treated = cf$tr,
                         n_post = cf$post, effect = 1, seed = 2026)
    Y <- matrix(sim$y, cf$N, cf$TT)
    W <- matrix(sim$d, cf$N, cf$TT)
    for (threads in c(1L, 4L)) {
      plan <- get("trop_prepare", ns)(
        Y, W, method = "local",
        grid = trop_grid(time = c(0, 0.5, 1), unit = c(0, 0.5, 1), nn = cf$nn),
        cv_cells = cf$cells, se = "bootstrap", n_boot = 4L, seed = 7,
        control = trop_control(threads = threads)
      )
      search <- function() {
        .Call(get("C_local_search", ns), plan$p$Y, plan$p$W, plan$p$obs, NULL,
              plan$grid$time, plan$grid$unit, plan$grid$nn, plan$cells,
              plan$max_cycles, plan$ctrl)
      }
      invisible(search())
      times <- double(cf$inner)
      for (j in seq_len(cf$inner)) {
        got <- search()
        times[j] <- got$kernel_seconds
      }
      est <- .Call(get("C_local_estimate", ns), plan$p$Y, plan$p$W, plan$p$obs,
                   NULL, got$lambda, plan$ctrl)
      boot <- .Call(get("C_bootstrap", ns), plan$p$Y, plan$p$W, plan$p$obs,
                    NULL, FALSE, 0L, 0L, got$lambda, plan$rows, plan$ctrl)
      results[[paste(cf$id, threads)]] <- list(
        config = cf$id, threads = threads, inner = cf$inner, seconds = mean(times),
        lambda = got$lambda, score = got$score, trace = got$trace,
        att = est$att, tau = est$tau, draws = boot$draws,
        n_skipped = est$n_skipped, n_nonconverged = est$n_nonconverged
      )
    }
  }
  saveRDS(results, output)
}

main <- function(args) {
  if (length(args) < 2L) stop("Supply the baseline and candidate library paths.")
  libs <- normalizePath(args[1:2], mustWork = TRUE)
  rounds <- if (length(args) >= 3L) as.integer(args[3L]) else 7L
  if (is.na(rounds) || rounds < 2L) stop("Use at least two rounds.")
  script <- normalizePath(file.path("inst", "bench", "time_cache_timings.R"))
  exe <- file.path(R.home("bin"), "Rscript")
  scratch <- tempfile(fileext = ".rds")
  on.exit(unlink(scratch), add = TRUE)
  raw <- list()
  max_diff <- function(a, b) {
    stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-8)))
    finite <- is.finite(a) & is.finite(b)
    if (any(finite)) max(abs(a[finite] - b[finite])) else 0
  }
  for (round in seq_len(rounds)) {
    pair <- list()
    # Alternate which revision runs first to limit order effects.
    for (revision in if (round %% 2L) c(1L, 2L) else c(2L, 1L)) {
      status <- system2(exe, c(shQuote(script), "--worker", shQuote(libs[revision]),
                              shQuote(scratch)))
      if (status != 0L) stop("A benchmark worker failed.")
      pair[[revision]] <- readRDS(scratch)
    }
    for (key in names(pair[[1L]])) {
      before <- pair[[1L]][[key]]
      after <- pair[[2L]][[key]]
      stopifnot(identical(before$lambda, after$lambda),
                identical(before$trace[, 1:3], after$trace[, 1:3]),
                identical(before$n_skipped, after$n_skipped),
                identical(before$n_nonconverged, after$n_nonconverged))
      score_diff <- max_diff(before$trace[, 4L], after$trace[, 4L])
      tau_diff <- max_diff(before$tau, after$tau)
      draw_diff <- max_diff(before$draws, after$draws)
      for (revision in 1:2) {
        x <- pair[[revision]][[key]]
        one <- pair[[revision]][[paste(x$config, 1L)]]
        stopifnot(identical(x$lambda, one$lambda), identical(x$trace, one$trace),
                  identical(x$tau, one$tau), identical(x$draws, one$draws))
        raw[[length(raw) + 1L]] <- data.frame(
          config = x$config, threads = x$threads, inner = x$inner, round = round,
          revision = c("baseline", "cache")[revision], seconds = x$seconds,
          lambda_time = x$lambda[1L], lambda_unit = x$lambda[2L], lambda_nn = x$lambda[3L],
          score = x$score, n_evals = nrow(x$trace), att = x$att,
          max_trace_score_diff = score_diff, max_tau_diff = tau_diff,
          max_bootstrap_diff = draw_diff, identical_across_threads = TRUE
        )
      }
    }
    utils::write.csv(do.call(rbind, raw),
                     file.path("inst", "bench", "time-cache-rounds.csv"), row.names = FALSE)
    message("Completed pair ", round, " of ", rounds)
  }
  raw <- do.call(rbind, raw)
  groups <- split(raw, interaction(raw$config, raw$threads))
  summary <- lapply(groups, function(x) {
    before <- x[x$revision == "baseline", ]
    after <- x[x$revision == "cache", ]
    lr <- log(before$seconds / after$seconds)
    margin <- stats::qt(0.975, length(lr) - 1L) * stats::sd(lr) / sqrt(length(lr))
    data.frame(config = x$config[1L], threads = x$threads[1L], rounds = rounds,
               baseline_median_s = stats::median(before$seconds),
               cache_median_s = stats::median(after$seconds),
               paired_speedup = exp(mean(lr)),
               speedup_lower_95 = exp(mean(lr) - margin),
               speedup_upper_95 = exp(mean(lr) + margin))
  })
  summary <- do.call(rbind, summary)
  utils::write.csv(summary, file.path("inst", "bench", "time-cache-summary.csv"),
                   row.names = FALSE)
  writeLines(c(trimws(capture.output(sessionInfo()), which = "right"),
               "Baseline: e85d03f (PR #4), recurrence",
               "Candidate: local mirrored direct-exponential time cache"),
             file.path("inst", "bench", "time-cache-session.txt"))
  print(summary, row.names = FALSE)
}

if (length(args) && args[1L] == "--worker") {
  worker(args[2L], args[3L])
} else {
  main(args)
}
