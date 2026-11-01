# Calculate the largest differences between fasttrop with covariates and the
# covariate references. The script data-raw/oracle_covariates.py made the
# references. Reference A is the global estimator of the authors' package with
# the term X phi added. Reference B is an independent solution of Eq. 13 for
# both methods. The script writes inst/validation/oracle-covariates-comparison.csv.
#
# Usage, in the package folder:
#   Rscript inst/validation/oracle_covariates_report.R
library(fasttrop)
fx_dir <- file.path("tests", "testthat", "fixtures")
read_fx <- function(nm) {
  d <- utils::read.csv(file.path(fx_dir, paste0(nm, ".csv")))
  N <- max(d$unit)
  TT <- max(d$time)
  idx <- cbind(d$unit, d$time)
  Y <- matrix(NA_real_, N, TT)
  W <- matrix(0L, N, TT)
  Y[idx] <- d$y
  W[idx] <- d$d
  X <- array(0, c(N, TT, 2), dimnames = list(NULL, NULL, c("x1", "x2")))
  X[cbind(idx, 1L)] <- d$x1
  X[cbind(idx, 2L)] <- d$x2
  list(data = d, Y = Y, W = W, X = X)
}
lam <- utils::read.csv(file.path(fx_dir, "lambdas.csv"))
lam_of <- function(cfg) {
  r <- lam[lam$config == cfg, ]
  c(time = r$time, unit = r$unit, nn = r$nn)
}
ctl <- trop_control(tol = 1e-10, max_iter = 100000L, threads = 4L)
fit_of <- function(fx, method, l, cm = "joint") {
  suppressWarnings(trop_matrix(fx$Y, fx$W, method = method, lambda = l, se = "none",
                               control = ctl, X = fx$X, covariate_method = cm))
}
rows <- list()
add <- function(quantity, method, reference, diff) {
  rows[[length(rows) + 1L]] <<- data.frame(quantity = quantity, method = method,
                                           reference = reference, n = length(diff),
                                           max_abs_diff = max(diff))
}

ref <- utils::read.csv(file.path(fx_dir, "oracle_cov.csv"))
ref$att_ft <- ref$phi1_ft <- ref$phi2_ft <- NA_real_
for (k in seq_len(nrow(ref))) {
  fit <- fit_of(read_fx(ref$panel[k]), ref$method[k], lam_of(ref$config[k]))
  ref$att_ft[k] <- fit$estimate
  ref$phi1_ft[k] <- fit$covariates$coef[["x1"]]
  ref$phi2_ft[k] <- if ("x2" %in% names(fit$covariates$coef)) fit$covariates$coef[["x2"]] else NA
}
for (m in unique(ref$method)) for (r in unique(ref$reference[ref$method == m])) {
  s <- ref[ref$method == m & ref$reference == r, ]
  add("ATT", m, r, abs(s$att_ft - s$att))
  if (m == "global") {
    add("phi", m, r, c(abs(s$phi1_ft - s$phi_1), abs(s$phi2_ft - s$phi_2)[!is.na(s$phi_2)]))
  }
}

cells <- utils::read.csv(file.path(fx_dir, "oracle_cov_cells.csv"))
cells <- cells[cells$method == "local", ]
d_tau <- double(0)
for (key in unique(paste(cells$panel, cells$config))) {
  r <- cells[paste(cells$panel, cells$config) == key, ]
  fx <- read_fx(r$panel[1])
  fit <- suppressWarnings(trop(fx$data, "y", "d", "unit", "time", method = "local",
                               lambda = lam_of(r$config[1]), se = "none", control = ctl,
                               covariates = c("x1", "x2")))
  m <- merge(r, as.data.frame(fit), by = c("unit", "time"), suffixes = c(".ref", ".fit"))
  d_tau <- c(d_tau, abs(m$tau.fit - m$tau.ref))
}
add("effect of each treated cell", "local", "B", d_tau)

res <- utils::read.csv(file.path(fx_dir, "oracle_cov_resid.csv"))
d_phi <- double(0)
for (k in seq_len(nrow(res))) {
  fit <- fit_of(read_fx(res$panel[k]), "local", c(time = 0, unit = 0, nn = Inf), "residualize")
  d_phi <- c(d_phi, abs(fit$covariates$coef[["x1"]] - res$phi_1[k]))
  if (!is.na(res$phi_2[k])) d_phi <- c(d_phi, abs(fit$covariates$coef[["x2"]] - res$phi_2[k]))
}
add("phi, residualize", "both", "TWFE regression", d_phi)

out <- do.call(rbind, rows)
utils::write.csv(out, file.path("inst", "validation", "oracle-covariates-comparison.csv"),
                 row.names = FALSE)
print(out, digits = 3)
