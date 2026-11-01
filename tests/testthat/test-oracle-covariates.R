# Agreement with the covariate references for fixed tuning parameters. The
# script data-raw/oracle_covariates.py made the references on the panels
# fixtures/XC1.csv, XC2.csv, and XC3.csv. fixtures/oracle_cov_meta.json gives
# the commits and the solver options.
#
# * Reference A is the global estimator of the authors' package with the term
#   X phi added to the model. It uses the conic solver CLARABEL. Thus, its
#   results must agree to 1e-4.
# * Reference B is an independent solution of Eq. 13 for the local and the
#   global methods. A proximal step after CLARABEL makes it accurate. Thus, its
#   results must agree to 1e-6.
# * In XC3, x2 does not change in time. The unit effects absorb it.

cov_fixture <- function(name) {
  d <- utils::read.csv(test_path("fixtures", paste0(name, ".csv")))
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

cov_fit <- function(fx, method, l, covariate_method = "joint") {
  suppressWarnings(trop_matrix(fx$Y, fx$W, method = method, lambda = l, se = "none",
                               control = tight_control(), X = fx$X,
                               covariate_method = covariate_method))
}

test_that("covariate ATTs and coefficients match references A and B", {
  ref <- utils::read.csv(test_path("fixtures", "oracle_cov.csv"))
  lam <- fixture_lambdas()
  for (k in seq_len(nrow(ref))) {
    fx <- cov_fixture(ref$panel[k])
    l <- lambda_of(lam[lam$config == ref$config[k], ])
    fit <- cov_fit(fx, ref$method[k], l)
    tol <- if (ref$reference[k] == "A") 1e-4 else 1e-6
    lab <- paste(ref$panel[k], ref$config[k], ref$method[k], "vs", ref$reference[k])
    expect_max_abs_diff(fit$estimate, ref$att[k], tol, label = paste(lab, "ATT"))
    if (ref$method[k] == "global") {
      expect_max_abs_diff(fit$covariates$coef[["x1"]], ref$phi_1[k], tol, label = paste(lab, "phi_1"))
      if (is.na(ref$phi_2[k])) {
        expect_false("x2" %in% fit$covariates$used)
      } else {
        expect_max_abs_diff(fit$covariates$coef[["x2"]], ref$phi_2[k], tol,
                            label = paste(lab, "phi_2"))
      }
    }
  }
})

test_that("local covariate cell effects match reference B cell by cell", {
  ref <- utils::read.csv(test_path("fixtures", "oracle_cov_cells.csv"))
  ref <- ref[ref$method == "local", ]
  lam <- fixture_lambdas()
  for (key in unique(paste(ref$panel, ref$config))) {
    r <- ref[paste(ref$panel, ref$config) == key, ]
    fx <- cov_fixture(r$panel[1])
    fit <- suppressWarnings(trop(fx$data, "y", "d", "unit", "time", method = "local",
                                 lambda = lambda_of(lam[lam$config == r$config[1], ]),
                                 se = "none", control = tight_control(),
                                 covariates = c("x1", "x2")))
    got <- as.data.frame(fit)
    m <- merge(r, got, by = c("unit", "time"), suffixes = c(".ref", ".fit"))
    expect_equal(nrow(m), nrow(r), label = paste(key, "matched cells"))
    expect_max_abs_diff(m$tau.fit, m$tau.ref, 1e-6, label = paste(key, "tau"))
    expect_max_abs_diff(fit$covariates$coef[["x1"]], mean(r$phi_1), 1e-6,
                        label = paste(key, "mean phi_1"))
  }
})

test_that("the residualize coefficients match the reference", {
  ref <- utils::read.csv(test_path("fixtures", "oracle_cov_resid.csv"))
  l <- c(time = 0, unit = 0, nn = Inf)
  for (k in seq_len(nrow(ref))) {
    fx <- cov_fixture(ref$panel[k])
    fit <- cov_fit(fx, "local", l, "residualize")
    expect_max_abs_diff(fit$covariates$coef[["x1"]], ref$phi_1[k], 1e-8,
                        label = paste(ref$panel[k], "phi_1"))
    if (!is.na(ref$phi_2[k])) {
      expect_max_abs_diff(fit$covariates$coef[["x2"]], ref$phi_2[k], 1e-8,
                          label = paste(ref$panel[k], "phi_2"))
    }
  }
})
