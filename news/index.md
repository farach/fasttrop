# Changelog

## fasttrop (development version)

- Time weights now use direct exponential values, with zero after
  underflow. This corrects the identifiability test for distant cells in
  long panels. The correction applies to local estimation,
  cross-validation, and bootstrap draws.

## fasttrop 0.2.0

### Covariates

- [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) has
  the new argument `covariates`, and
  [`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
  has the new argument `X`. They add the term (X\_{it}’) of Eq. 13 in
  Section 6.2 of the paper.
- The new argument `covariate_method` selects how to calculate (). With
  `"joint"` (the default), each weighted fit calculates () together with
  the other parts of the model. With `"residualize"`, one TWFE
  regression on the untreated cells calculates () first.
- The global method can also use covariates. This is an extension of
  fasttrop. Section 6.2 of the paper adds covariates to the local
  method. The Python package from the authors does not have covariates.
- The result has the new element `covariates` with the coefficients, the
  method, and the names of the covariates that the fit used.
- The fixed effects absorb a covariate that does not change in time or
  does not change between units. The function gives a warning message
  for such a covariate.

### Validation

- New tests compare the results with covariates with two references. The
  script `data-raw/oracle_covariates.py` adds the covariate term to the
  code of the authors and also solves the problem with a second method.
- The new script `inst/validation/covariates.R` runs a simulation with a
  known effect and a covariate that is different for the treated units.
- The new scripts in `inst/validation/stata` compare fasttrop with the
  Stata package `trop`. They only call that package.
- The new script `inst/bench/covariate_timings.R` measures the time cost
  of covariates.
- The new article “Use covariates” explains the two procedures and the
  results.

## fasttrop 0.1.0

This is the first version of the package.

### Estimator

- The package calculates the Triply Robust Panel (TROP) estimator of
  Athey, Imbens, Qu, and Viviano (2026).
- The local method (`method = "local"`) is the estimator of the paper
  for each treated cell. The global method (`method = "global"`) is the
  block estimator of the Python package from the authors.
- Leave-one-out cross-validation or placebo cross-validation selects the
  tuning parameters.
- The unit bootstrap of Algorithm 3 gives standard errors. The default
  confidence interval is a normal interval. A percentile interval is
  also available.

### Speed

- The model fits run in C with OpenMP and with the BLAS and LAPACK
  libraries of R. The number of threads does not change the results.
- On Windows, the C code tells Windows not to use power throttling on
  its threads while it runs. Power throttling can make a long
  calculation in an R session slower. To keep the usual Windows policy,
  use `trop_control(throttle = TRUE)`.
- The function
  [`trop_threads()`](https://farach.github.io/fasttrop/reference/trop_threads.md)
  gives at most 2 threads when `R CMD check` sets the environment
  variable `_R_CHECK_LIMIT_CORES_`.

### Documentation

- The package has a website with articles in plain words:
  <https://farach.github.io/fasttrop/>.
- The vignette “Get started with fasttrop” shows a first analysis step
  by step.
- `citation("fasttrop")` gives the references for the paper and the
  package.

### Validation

- The documents use the name “TWFE imputation” for the special case
  `lambda = c(time = 0, unit = 0, nn = Inf)`. This method fits a two-way
  fixed effects model to the untreated cells only. In a staggered
  design, it is different from a conventional TWFE regression.
- Two new tests compare TWFE imputation with a conventional TWFE
  regression. The test suite has 478 expectations.
- The new script `inst/validation/staggered.R` runs simulations with a
  staggered design and a known effect. It compares a conventional TWFE
  regression, TWFE imputation, and the local method of TROP. The folder
  `inst/validation` has the saved results.
- The article “How we checked fasttrop” has the results of the new
  simulations. The article “What TROP can and cannot do” has the limits.
