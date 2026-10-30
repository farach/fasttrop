# Report the largest differences between fasttrop and each reference (C1).
library(fasttrop)
fx_dir <- file.path("tests", "testthat", "fixtures")
read_fx <- function(nm) {
  df <- utils::read.csv(file.path(fx_dir, paste0(nm, ".csv")))
  pm <- fasttrop:::panel_matrices(df, "y", "d", "unit", "time")
  list(data = df, Y = pm$Y, W = pm$W)
}
lam <- utils::read.csv(file.path(fx_dir, "lambdas.csv"))
lam_of <- function(cfg) { r <- lam[lam$config == cfg, ]; c(time = r$time, unit = r$unit, nn = r$nn) }
ctl <- trop_control(tol = 1e-10, max_iter = 100000L, threads = 4L)

att <- utils::read.csv(file.path(fx_dir, "oracle.csv"))
att$fasttrop <- NA_real_
t0 <- Sys.time()
for (k in seq_len(nrow(att))) {
  fx <- read_fx(att$dataset[k])
  att$fasttrop[k] <- suppressWarnings(trop_matrix(fx$Y, fx$W, method = att$method[k],
    lambda = lam_of(att$config[k]), se = "none", control = ctl))$estimate
}
att$diff <- abs(att$fasttrop - att$att)
cells <- utils::read.csv(file.path(fx_dir, "oracle_cells.csv"))
cells$fasttrop <- NA_real_
for (key in unique(paste(cells$dataset, cells$config, cells$method))) {
  idx <- which(paste(cells$dataset, cells$config, cells$method) == key)
  r <- cells[idx[1], ]
  fx <- read_fx(r$dataset)
  fit <- suppressWarnings(trop(fx$data, "y", "d", "unit", "time", method = r$method,
                               lambda = lam_of(r$config), se = "none", control = ctl))
  got <- as.data.frame(fit)
  cells$fasttrop[idx] <- got$tau[match(paste(cells$unit[idx], cells$time[idx]),
                                       paste(got$unit, got$time))]
}
cells$diff <- abs(cells$fasttrop - cells$tau)

obj <- utils::read.csv(file.path(fx_dir, "oracle_objective.csv"))
sol <- utils::read.csv(file.path(fx_dir, "oracle_solver.csv"))
obj$fasttrop <- obj$rel_diff <- obj$max_fitted_diff <- NA_real_
for (k in seq_len(nrow(obj))) {
  o <- obj[k, ]
  s <- sol[sol$dataset == o$dataset & sol$config == o$config & sol$method == o$method, ]
  fx <- read_fx(o$dataset)
  idx <- cbind(match(s$unit, as.integer(rownames(fx$Y))), match(s$time, as.integer(colnames(fx$Y))))
  w <- matrix(0, nrow(fx$Y), ncol(fx$Y)); w[idx] <- s$w
  rf <- matrix(NA_real_, nrow(fx$Y), ncol(fx$Y)); rf[idx] <- s$fitted
  Y <- unname(fx$Y)
  fit <- fasttrop:::trop_fit_weighted(Y, w, o$lambda_nn, ctl)
  fitted <- outer(fit$alpha, rep(1, ncol(Y))) + outer(rep(1, nrow(Y)), fit$beta) + fit$L
  obj$fasttrop[k] <- sum(w * (Y - fitted)^2) + o$lambda_nn * sum(svd(fit$L)$d)
  obj$max_fitted_diff[k] <- max(abs(fitted[w > 0] - rf[w > 0]))
}
obj$rel_diff <- (obj$fasttrop - obj$objective) / obj$objective
cat(sprintf("elapsed %.1f s\n", as.double(difftime(Sys.time(), t0, units = "secs"))))
cat("\nATT: max |fasttrop - reference| by method and reference\n")
print(aggregate(diff ~ method + reference, data = att, FUN = max), digits = 3)
cat("\nATT differences by configuration (max over datasets)\n")
print(aggregate(diff ~ method + reference + config, data = att, FUN = max), digits = 3)
cat("\nCell effects vs diff-diff: max |difference| by method\n")
print(aggregate(diff ~ method, data = cells, FUN = max), digits = 3)
cat("\nSolver on identical weights: objective (fasttrop - diff-diff) / diff-diff, fitted values\n")
print(obj[, c("dataset", "config", "method", "objective", "fasttrop", "rel_diff", "max_fitted_diff")],
      digits = 6)
utils::write.csv(att, file.path("inst", "validation", "oracle-att-comparison.csv"), row.names = FALSE)
utils::write.csv(obj, file.path("inst", "validation", "oracle-objective-comparison.csv"), row.names = FALSE)
