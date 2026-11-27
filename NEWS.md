# fasttrop (development version)

## Effects by group

- `trop_effects()` accepts a named list in `by`. It makes all the sets of
  groups from one set of bootstrap draws and returns a named list of
  `trop_effects` objects. Give this list to `trop_stepdown()`. Before, each
  set of groups needed a new set of draws.
- The new argument `anticipation` of `trop_effects()` gives event time 0 to
  the first period after the anticipation periods. The cohort is this period.
  Use it when you set the anticipation periods to `W = 1` and remove them from
  `target`. The default, 0, gives the result of version 0.3.0.
- The bootstrap of `trop_effects()` runs in the C code, in parallel with
  `threads`. Before, it ran in R, one draw after another.

## Speed

- The local estimator calculates time-decay weights once for each fit or
  cross-validation call and shares them across worker threads.
  It still uses direct exponential values to preserve zero weights after
  underflow.

## Bug fixes

- In `trop_effects()`, a treated cell that `target` removes is not in a group.
  Before, such cells made a group with no cells. That group had no
  estimate, and `trop_stepdown()` then gave an error.

# fasttrop 0.3.0

This version adds tools for real studies with many treated cells, many groups,
and many tests. The default of each new argument gives the same result as
version 0.2.0.

## Effects by group

- The new function `trop_effects()` gives the mean effect of groups of treated
  cells. The argument `by` selects the groups: `"event_time"`, `"cohort"`,
  `"time"`, `"unit"`, `"overall"`, or a matrix with one group number for each
  cell. It returns the estimates, the standard errors, the intervals, and a
  simultaneous (sup-t) band for all groups.
- The new methods `print()` and `as.data.frame()` work for the result of
  `trop_effects()`.
- `trop_effects()` uses the bootstrap rows of the fit. Thus, the draws of all
  groups come from the same resamples. For the group `"overall"`, the draws are
  the same as the draws of the fit.

## Many tests

- The new function `trop_boot_rows()` makes the bootstrap rows once. Give the
  rows to the argument `boot_rows` of `trop()` and `trop_matrix()` to use the
  same resamples for many fits.
- The new function `trop_stepdown()` adjusts p-values with the stepdown method
  of Romano and Wolf (2005). It accepts fits or the results of
  `trop_effects()`.
- The new argument `strata` resamples the units inside groups. The new argument
  `resample = "controls"` of `trop_boot_rows()` resamples only the untreated
  units.

## Target cells

- The new argument `target` of `trop()` and `trop_matrix()` keeps treated cells
  out of the estimate. A cell outside of the target stays in the panel, but the
  estimate does not use it. The fit counts these cells in
  `n_excluded_cells`. A target works with the local method only.

## Tuning

- The new function `trop_grid_scaled()` puts the tuning grid on the scale of
  the outcome. The setting `grid = "auto"` in `trop()` and `trop_matrix()` uses
  it with the default grid.
- Cross-validation now gives a warning message when it chooses a value at the
  edge of the grid. The result `cv$boundary` names these values.
- The new argument `demean = TRUE` subtracts the mean of each unit before the
  fit.

## Records

- A fit now stores the version of fasttrop, and the commit for an installation
  from GitHub, in `fit$package`. `summary()` prints them.
- A fit now stores the panel in `fit$panel` and the bootstrap rows in
  `fit$boot$rows`. A fit from an older version does not have them, so
  `trop_effects()` cannot use such a fit.
- `print()` shows the number of treated cells that are outside of the target.

## Fixes

- Time weights now use direct exponential values, with zero after underflow.
  This corrects the identifiability test for distant cells in long panels.
  The correction applies to local estimation, cross-validation, and bootstrap
  draws.

# fasttrop 0.2.0

## Covariates

- `trop()` has the new argument `covariates`, and `trop_matrix()` has the new
  argument `X`. They add the term \(X_{it}'\phi\) of Eq. 13 in Section 6.2 of
  the paper.
- The new argument `covariate_method` selects how to calculate \(\phi\). With
  `"joint"` (the default), each weighted fit calculates \(\phi\) together with
  the other parts of the model. With `"residualize"`, one TWFE regression on
  the untreated cells calculates \(\phi\) first.
- The global method can also use covariates. This is an extension of fasttrop.
  Section 6.2 of the paper adds covariates to the local method. The Python
  package from the authors does not have covariates.
- The result has the new element `covariates` with the coefficients, the
  method, and the names of the covariates that the fit used.
- The fixed effects absorb a covariate that does not change in time or does
  not change between units. The function gives a warning message for such a
  covariate.

## Validation

- New tests compare the results with covariates with two references. The
  script `data-raw/oracle_covariates.py` adds the covariate term to the code of
  the authors and also solves the problem with a second method.
- The new script `inst/validation/covariates.R` runs a simulation with a known
  effect and a covariate that is different for the treated units.
- The new scripts in `inst/validation/stata` compare fasttrop with the Stata
  package `trop`. They only call that package.
- The new script `inst/bench/covariate_timings.R` measures the time cost of
  covariates.
- The new article "Use covariates" explains the two procedures and the
  results.

# fasttrop 0.1.0

This is the first version of the package.

## Estimator

- The package calculates the Triply Robust Panel (TROP) estimator of Athey,
  Imbens, Qu, and Viviano (2026).
- The local method (`method = "local"`) is the estimator of the paper for each
  treated cell. The global method (`method = "global"`) is the block estimator
  of the Python package from the authors.
- Leave-one-out cross-validation or placebo cross-validation selects the
  tuning parameters.
- The unit bootstrap of Algorithm 3 gives standard errors. The default
  confidence interval is a normal interval. A percentile interval is also
  available.

## Speed

- The model fits run in C with OpenMP and with the BLAS and LAPACK libraries
  of R. The number of threads does not change the results.
- On Windows, the C code tells Windows not to use power throttling on its
  threads while it runs. Power throttling can make a long calculation in an R
  session slower. To keep the usual Windows policy, use
  `trop_control(throttle = TRUE)`.
- The function `trop_threads()` gives at most 2 threads when `R CMD check`
  sets the environment variable `_R_CHECK_LIMIT_CORES_`.

## Documentation

- The package has a website with articles in plain words:
  <https://farach.github.io/fasttrop/>.
- The vignette "Get started with fasttrop" shows a first analysis step by
  step.
- `citation("fasttrop")` gives the references for the paper and the package.

## Validation

- The documents use the name "TWFE imputation" for the special case
  `lambda = c(time = 0, unit = 0, nn = Inf)`. This method fits a two-way fixed
  effects model to the untreated cells only. In a staggered design, it is
  different from a conventional TWFE regression.
- Two new tests compare TWFE imputation with a conventional TWFE regression.
  The test suite has 478 expectations.
- The new script `inst/validation/staggered.R` runs simulations with a
  staggered design and a known effect. It compares a conventional TWFE
  regression, TWFE imputation, and the local method of TROP. The folder
  `inst/validation` has the saved results.
- The article "How we checked fasttrop" has the results of the new
  simulations. The article "What TROP can and cannot do" has the limits.
