# Generates the oracle fixture panels used by tests/testthat/test-oracle.R.
# Run from the package root: Rscript data-raw/make_fixtures.R
# Then run data-raw/oracle_fixtures.py to compute the reference answers.
devtools::load_all(".", quiet = TRUE)

out <- file.path("tests", "testthat", "fixtures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

specs <- list(
  D1 = list(n_units = 12, n_periods = 10, n_treated = 3, n_post = 3, design = "block", seed = 101),
  D2 = list(n_units = 20, n_periods = 15, n_treated = 4, n_post = 4, design = "block", seed = 202),
  D3 = list(n_units = 15, n_periods = 12, n_treated = 5, n_post = 4, design = "staggered", seed = 303),
  D4 = list(n_units = 16, n_periods = 12, n_treated = 1, n_post = 3, design = "block", seed = 404)
)

for (nm in names(specs)) {
  s <- specs[[nm]]
  d <- trop_simulate(n_units = s$n_units, n_periods = s$n_periods, n_treated = s$n_treated,
                     n_post = s$n_post, effect = 1, design = s$design, seed = s$seed)
  utils::write.csv(d[, c("unit", "time", "y", "d")], file.path(out, paste0(nm, ".csv")),
                   row.names = FALSE)
}

lambdas <- data.frame(
  config = c("twfe", "weighted_twfe", "mc", "trop_a", "trop_b"),
  time = c(0, 0.5, 0, 0.3, 1.0),
  unit = c(0, 0.5, 0, 0.8, 0.2),
  nn = c(Inf, Inf, 0.1, 0.05, 1.0)
)
utils::write.csv(lambdas, file.path(out, "lambdas.csv"), row.names = FALSE)
cat("Wrote", length(specs), "panels and", nrow(lambdas), "lambda configurations to", out, "\n")
