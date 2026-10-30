# Matched timings against the reference implementations (council condition C4,
# open defeater D-A-O1).
#
# Same panels, same fixed tuning parameters, and each package's default solver
# settings. fasttrop runs with 1 and 4 threads; cfcompare (pure R) with its
# per-cell estimator; the Python references are timed by reference_timings.py,
# which this script runs first if the Python results are missing.
#
# Usage (from the package root):
#   Rscript inst/bench/reference_timings.R            # fasttrop and cfcompare
#   python inst/bench/reference_timings.py <trop> <diff-diff>
#   Rscript inst/bench/reference_timings.R combine    # merge into one table
library(fasttrop)
bench <- file.path("inst", "bench")
lam <- c(time = 0.3, unit = 0.8, nn = 0.5)
reps <- 3L
args <- commandArgs(trailingOnly = TRUE)

make_panel <- function(N, TT, ntr, post, seed) {
  sim <- trop_simulate(n_units = N, n_periods = TT, n_treated = ntr, n_post = post,
                       effect = 1, seed = seed)
  sim[, c("unit", "time", "y", "d")]
}
panels <- list(panel_20x15 = make_panel(20, 15, 4, 4, 11), panel_40x20 = make_panel(40, 20, 8, 5, 12))

if (!identical(args, "combine")) {
  for (nm in names(panels)) {
    utils::write.csv(panels[[nm]], file.path(bench, paste0(nm, ".csv")), row.names = FALSE)
  }
  best_of <- function(f, n = reps) {
    out <- NULL
    tm <- min(vapply(seq_len(n), function(i) {
      t0 <- bench::hires_time()
      out <<- f()
      as.double(bench::hires_time() - t0)
    }, double(1)))
    list(seconds = tm, value = out)
  }
  rows <- list()
  add <- function(panel, task, impl, r) {
    rows[[length(rows) + 1L]] <<- data.frame(panel = panel, task = task, implementation = impl,
                                             seconds = r$seconds, value = r$value)
  }
  for (nm in names(panels)) {
    d <- panels[[nm]]
    pm <- fasttrop:::panel_matrices(d, "y", "d", "unit", "time")
    for (th in c(1L, 4L)) {
      ctl <- trop_control(threads = th)
      impl <- sprintf("fasttrop (C, %d thread%s)", th, if (th > 1L) "s" else "")
      add(nm, "local ATT, fixed lambda", impl, best_of(function() {
        trop_matrix(pm$Y, pm$W, method = "local", lambda = lam, se = "none", control = ctl)$estimate
      }))
      add(nm, "global ATT, fixed lambda", impl, best_of(function() {
        trop_matrix(pm$Y, pm$W, method = "global", lambda = lam, se = "none", control = ctl)$estimate
      }))
      if (nm == "panel_20x15") {
        inp <- fasttrop:::as_panel_inputs(pm$Y, pm$W)
        add(nm, "LOOCV Q(lambda), all untreated cells", impl, best_of(function() {
          .Call(fasttrop:::C_local_loocv, inp$Y, inp$W, inp$obs, lam, NULL, unclass(ctl))$Q
        }))
      }
    }
    if (requireNamespace("cfcompare", quietly = TRUE)) {
      add(nm, "local ATT, fixed lambda", "cfcompare (pure R, per_cell)", best_of(function() {
        cfcompare::trop(d, "y", "d", "unit", "time", lambda = as.list(lam), anchor = "per_cell",
                        se = "none")$estimate
      }, n = 1L))
    }
    message(nm, " done")
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(bench, "reference-timings-r.csv"), row.names = FALSE)
  print(res, digits = 6)
} else {
  r <- utils::read.csv(file.path(bench, "reference-timings-r.csv"))
  py <- utils::read.csv(file.path(bench, "reference-timings-python.csv"))
  all <- rbind(r, py)
  base <- all[all$implementation == "fasttrop (C, 1 thread)", c("panel", "task", "seconds")]
  names(base)[3] <- "fasttrop_1thread_s"
  all <- merge(all, base, by = c("panel", "task"))
  all$slower_than_fasttrop_1thread <- all$seconds / all$fasttrop_1thread_s
  all <- all[order(all$panel, all$task, all$seconds), ]
  utils::write.csv(all, file.path(bench, "reference-timings.csv"), row.names = FALSE)
  print(all[, c("panel", "task", "implementation", "seconds", "value",
                "slower_than_fasttrop_1thread")], digits = 4, row.names = FALSE)
}
