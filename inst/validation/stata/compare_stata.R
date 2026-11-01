# Compare the fasttrop covariate results with the results of the Stata package
# gorgeousfish/TROP. First run trop_covariates.do in Stata. Then run this
# script in the fasttrop source folder:
#   Rscript inst/validation/stata/compare_stata.R
#
# The script writes inst/validation/stata/stata-comparison.csv. For each
# panel, configuration, and method, it gives the ATT of the two packages
# without covariates and with covariates, and the covariate coefficients.
# The column "shift_diff" is the difference between the two packages in the
# change of the ATT that the covariates make. If the two packages use
# different conventions for the estimator without covariates, this column
# shows only the difference that comes from the covariate step.
#
# The two packages use different numerical methods. Thus, small differences
# are possible. The Stata package gives one coefficient vector for the local
# method. fasttrop gives the mean of the per-cell coefficients. Thus, compare
# the local coefficients with care.
library(fasttrop)
dir <- file.path("inst", "validation", "stata")
fix <- file.path("tests", "testthat", "fixtures")
st <- utils::read.csv(file.path(dir, "stata-covariates.csv"))
lam <- utils::read.csv(file.path(fix, "lambdas.csv"))
ctl <- trop_control(tol = 1e-10, max_iter = 100000L)

fit_one <- function(panel, config, method, covariates) {
  d <- utils::read.csv(file.path(fix, paste0(panel, ".csv")))
  l <- lam[lam$config == config, ]
  cv <- if (covariates) c("x1", "x2") else NULL
  f <- tryCatch(suppressWarnings(trop(d, "y", "d", "unit", "time", method = method,
                                      lambda = c(time = l$time, unit = l$unit, nn = l$nn),
                                      se = "none", control = ctl, covariates = cv)),
                error = function(e) NULL)
  if (is.null(f)) return(c(att = NA, g1 = NA, g2 = NA))
  co <- if (covariates) f$covariates$coef else c(x1 = NA, x2 = NA)
  c(att = f$estimate, g1 = unname(co["x1"]), g2 = unname(co["x2"]))
}

ft <- t(mapply(fit_one, st$panel, st$config, st$method, st$covariates == 1))
st$att_fasttrop <- ft[, "att"]
st$gamma1_fasttrop <- ft[, "g1"]
st$gamma2_fasttrop <- ft[, "g2"]
st$att_diff <- st$att_fasttrop - st$att

key <- paste(st$panel, st$config, st$method)
base <- st[st$covariates == 0, ]
cov <- st[st$covariates == 1, ]
cov <- cov[match(paste(base$panel, base$config, base$method), key[st$covariates == 1]), ]
out <- data.frame(panel = base$panel, config = base$config, method = base$method,
                  att_stata = base$att, att_fasttrop = base$att_fasttrop,
                  att_cov_stata = cov$att, att_cov_fasttrop = cov$att_fasttrop,
                  shift_diff = (cov$att_fasttrop - base$att_fasttrop) - (cov$att - base$att),
                  gamma1_stata = cov$gamma1, gamma1_fasttrop = cov$gamma1_fasttrop,
                  gamma2_stata = cov$gamma2, gamma2_fasttrop = cov$gamma2_fasttrop,
                  stata_rc = cov$rc)
utils::write.csv(out, file.path(dir, "stata-comparison.csv"), row.names = FALSE)
print(out, digits = 4)
ok <- is.finite(out$att_cov_stata)
message(sprintf("Largest |ATT difference| without covariates: %.2g",
                max(abs(out$att_fasttrop - out$att_stata)[ok], na.rm = TRUE)))
message(sprintf("Largest |ATT difference| with covariates: %.2g",
                max(abs(out$att_cov_fasttrop - out$att_cov_stata)[ok], na.rm = TRUE)))
message(sprintf("Largest |shift difference|: %.2g", max(abs(out$shift_diff)[ok], na.rm = TRUE)))
