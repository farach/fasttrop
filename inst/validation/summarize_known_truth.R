# Make the summaries of the simulations from the saved draws. Thus, all
# simulations give the same statistics. Run this script in the package folder
# after known_truth.R.
tau <- 1
z <- stats::qnorm(0.975)
dir <- file.path("inst", "validation")
files <- list.files(dir, pattern = "^known-truth-draws-.*\\.csv$", full.names = TRUE)
for (f in files) {
  d <- utils::read.csv(f)
  summ <- do.call(rbind, lapply(split(d, d$method), function(x) {
    cov_n <- mean(abs(x$estimate - tau) <= z * x$se)
    cov_p <- mean(x$lo <= tau & tau <= x$hi)
    data.frame(design = x$design[1], method = x$method[1], reps = nrow(x),
               bias = mean(x$estimate) - tau, mc_se_bias = stats::sd(x$estimate) / sqrt(nrow(x)),
               sd = stats::sd(x$estimate), mean_se = mean(x$se),
               rmse = sqrt(mean((x$estimate - tau)^2)),
               coverage_normal = cov_n, mc_se_normal = sqrt(cov_n * (1 - cov_n) / nrow(x)),
               coverage_percentile = cov_p,
               mc_se_percentile = sqrt(cov_p * (1 - cov_p) / nrow(x)),
               median_seconds = stats::median(x$seconds))
  }))
  summ <- summ[match(c("twfe", "local", "global"), summ$method), ]
  row.names(summ) <- NULL
  out <- sub("draws", "summary", f)
  utils::write.csv(summ, out, row.names = FALSE)
  print(summ, digits = 3)
}
