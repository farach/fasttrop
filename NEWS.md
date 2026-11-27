# fasttrop 0.4.0

This version adds tools for long studies. It adds diagnostics for each
bootstrap draw, checkpoints, and weights. It also counts the cells of the
estimate before the bootstrap.

With the default of each new argument, most results are the same as in
version 0.3.0. These results can change:

- The cross-validation can select different tuning parameters when some
  tuning values do not give a score for all cells or all placebo cells.
- The global placebo cross-validation can select different tuning parameters
  when `placebo_size` is more than 1. This also occurs when a placebo unit has
  no outcome in the held-out periods.
- The global estimator does not estimate a cell that the untreated cells do
  not identify. A bootstrap draw of the global estimator uses the same
  treated units as the fit.
- The sup-t band of `trop_effects()` uses only the draws with a value for all
  groups.
- In a `trop()` fit, a period without a row does not stop a treatment run.
  Thus, the event times and the cohorts of `trop_effects()` can change.

Refer to the sections Changes and Bug fixes.

## Plan the estimand

- The new function `trop_rollout()` makes the treatment, the target, and the
  event time from the start date of each unit. With `anticipation`, it marks
  the periods before the start as treated. These periods are not controls, and
  by default they are not in the target.
- The new argument `event_time` of `trop()` and `trop_matrix()` gives the
  event time of each treated cell. The fit keeps these event times. Then
  `trop_effects()` and `trop_cells()` use them for the event times, the
  cohorts, and `window`. Use the column `event_time` of `trop_rollout()`. This
  is necessary when the panel has no row for the start of a run, or when the
  anticipation periods start before the panel.
- In a run of treated periods, the event time must increase by 1 from one
  period to the next. Otherwise, the function stops with an error.
- The new function `trop_cells()` counts the cells of each group of
  `trop_effects()`. It counts the cells in the estimate, the cells outside of
  the target, and the cells with no effect. It does not show the effects. Thus,
  you can examine the estimand before you see the results or calculate the
  bootstrap. If no treated cell is in a group, the result has no rows, and the
  print method shows where the treated cells are.
- The new arguments `weights` and `unit_weights` of `trop_effects()` set the
  weight of each cell in a group. With `weights = "unit"`, each unit has the
  same total weight in a group. With `unit_weights`, each unit has its own
  weight. The bootstrap uses the same weights.
- The weights are relative. In each group, the functions divide the weights
  by the largest weight of the units of that group. Thus, the estimate of a
  group does not depend on the weights of other units. In a group, the
  largest positive weight must be less than 1e200 times the smallest positive
  weight. Named weights must agree exactly with the unit names of the fit.

## Effects by group

- `trop_effects()` accepts a named list in `by`. It makes all the sets of
  groups from one set of bootstrap draws and returns a named list of
  `trop_effects` objects. Give this list to `trop_stepdown()`. Before, each
  set of groups needed a new set of draws.
- The new argument `anticipation` of `trop_effects()` gives event time 0 to
  the first period after the anticipation periods. The cohort is this period.
  Use it when you set the anticipation periods to `W = 1` and remove them from
  `target`. The default, 0, gives the result of version 0.3.0. If the fit has
  event times from `event_time`, `anticipation` must be 0. The function gives
  a warning when a run starts in the first period, because its anticipation
  periods can be before the panel.
- The bootstrap of `trop_effects()` runs in the C code, in parallel with
  `threads`. Before, it ran in R, one draw after another.
- The result of `trop_effects()` has the new columns `n_units` and `n_draws`.
  The column `n_draws` is the number of draws with a value for each group. The
  function gives a warning for a group with fewer than 2 draws with a value.
- The new methods `coef()` and `vcov()` work for `trop_effects` objects. The
  new function `trop_contrast()` calculates linear combinations of group
  estimates from the same draws, for example the mean of event times 0 to 3.
  Each combination uses only the estimates with a coefficient that is not 0.
  `trop_stepdown()` accepts its result.
- The new function `trop_honest_inputs()` makes the inputs of the 'HonestDiD'
  package from event-time estimates with a reference period. Its element
  `l_vec` gives the weights of the event times that the sensitivity analysis
  examines. By default, this is event time 0. Give `l_vec` to 'HonestDiD',
  because the default of 'HonestDiD' is the first event time after the
  reference. The function does not accept event times with gaps.

## Diagnostics

- Each fit now has `n_failed`, the number of fits that failed. A fit fails if
  the singular value decomposition fails or if the result is not finite. Such
  a cell has no effect, and the function gives a warning message. Before, a
  failed decomposition gave a low-rank term of 0 without a message.
- If `dgesdd` fails, the C code now uses the LAPACK routine `dgesvd` for the
  decomposition.
- The bootstrap now counts the fits and the cells of each draw. A fit has
  these counts in `boot$diagnostics`, and a `trop_effects` object has them in
  the attribute `diagnostics`. They give the fits that did not converge, the
  fits that failed, and the cells that the draw cannot identify. The column
  `n_dropped` counts the cells that have an effect in the fit, but that the
  draw cannot identify. The functions give a warning message when a draw has
  such fits or cells. `summary()` and `print()` show them.
- A draw with no treated cell is not a failed fit. The messages now tell the
  two cases apart. They also tell when a group has no standard error.
- The trace of the cross-validation in `cv$trace` now has three counts for
  each grid point. They are the cells with a score (`n_valid`), the fits that
  did not converge, and the fits that failed. The function gives a warning
  message when fits of the cross-validation do not converge or fail.
- `fit$timing` now has the time in the C code of each step: the tuning, the
  estimate, and the bootstrap. `summary()` shows them.

## Long runs

- The new arguments `progress` and `checkpoint` of `trop_effects()` show the
  progress of the bootstrap and save the draws to a file after each block. If
  the calculation stops, the same call continues from the file. The results
  are equal to the results of one run without a stop. A fingerprint makes sure
  that the file is for the same calculation.
- Before the bootstrap starts, the function makes sure that it can write to
  the folder of the file. The test and the save use new temporary files with
  short names. Thus, they do not change other files in the folder.
- The file has a hash of the results. The function stops with an error if
  the results in the file do not agree with the hash. It also stops if a
  result cannot be the result of a bootstrap draw.

## Changes

- The cross-validation now compares two tuning values first by the number of
  cells with a score. For the global method, these are the placebo cells with
  an effect. With the same number, the smaller score is better. Before, a
  tuning value with a score for fewer cells sometimes won, because its score
  used fewer cells. The selection changes only when the tuning values have
  scores for different numbers of cells.
- In `cv$trace`, the column `n_valid` gives this number of cells. The new
  element `cv$n_cells` gives the number of all such cells.
- If no tuning value has a usable score, the cross-validation now stops with
  an error. Before, it selected the first value of the grid. If the chosen
  values have a score for fewer cells than the criterion has, the function
  gives a warning.
- The sup-t band of `trop_effects()` now uses only the draws with a value for
  all groups, as `trop_stepdown()` does. The centers and the standard errors
  of the band also come from these draws. Before, a draw with no value for
  some groups used the largest value of the other groups. Thus, the band was
  sometimes too narrow. When all draws have a value for all groups, the band
  does not change.
- The attribute `n_complete` gives the number of draws of the band. If fewer
  than 2 draws have a value for all groups, the band is not available, and the
  print method gives the reason.
- In `trop()`, a cell without a row now gets a treatment from the rows of its
  unit. It gets 1 if the nearest row before it has 1, and the nearest row after
  it has 1 or does not exist. Otherwise, it gets 0. Thus, a period without a
  row does not stop a treatment run.
- Before, a period without a row in a treatment run gave the design
  "general", a warning that the treatment stops, and an error with
  `method = "global"`. Also, `trop_effects()` started the event times and the
  cohort again after that period. The estimates do not change, because a cell
  without a row has no outcome.
- With `event_time`, the treated rows after a cell without a row can show
  that their run starts after the cell. Then the cell gets `W = 0`, and a new
  run starts after it. A run that starts before the cell must agree with the
  rows before the cell. Otherwise, the function stops with an error.
- If the panel has no rows just before the first treated row of a run, the
  start of the run is not certain. Then `trop()` gives a warning. Give
  `event_time` to get the correct event times.
- `trop_stepdown()` now returns an object of class `trop_stepdown`. Its print
  method shows the number of draws that it uses. The function gives a warning
  message when some draws do not have a value for all estimates.
- `summary()` now shows the number of draws without an estimate. Before, it
  called these draws failed draws. A failed fit is a different count.
- The functions now give an error when the rows of a `trop_effects` object do
  not match its draws. This occurs when you remove rows of the object.

## Speed and memory

- The local estimator now keeps the distance sums only for the units that are
  targets of a fit. Before, it kept the sums for all pairs of units, in each
  thread. For a panel with many untreated units, this uses much less memory.
  With `unit = 0`, the estimator keeps no distance sums. The estimates and the
  draws are the same as before.
- A bootstrap draw of the local estimator now uses the memory of its thread
  again. Before, each draw allocated new memory.
- The local estimator calculates time-decay weights once for each fit or
  cross-validation call and shares them across worker threads.
  It still uses direct exponential values to preserve zero weights after
  underflow.

## Bug fixes

- In `trop_effects()`, a treated cell that `target` removes is not in a group.
  Before, such cells made a group with no cells. That group had no
  estimate, and `trop_stepdown()` then gave an error.
- The global estimator now finds the treated cells that the untreated cells
  do not identify. An example is a period with no untreated outcome. It does
  not estimate these cells, and it counts them in `n_skipped`. Before, it gave
  a value for such cells without a message. This also applies to each
  bootstrap draw.
- `trop_stepdown()` and `trop_contrast()` now make sure that draw b of each
  object uses the same units. If each object has a different name for each
  unit, they compare the units. Before, they compared only the row numbers.
  Thus, two fits with different orders of the units gave wrong joint p-values
  without a message.
- If only some objects have valid unit names, draw b of all objects must have
  the same row numbers, and the functions give a warning. A valid name is not
  `NA`, and each unit has a different name. If no object has unit names, the
  functions compare the row numbers without a warning. The order of the
  objects does not change the result.
- A bootstrap draw of the global estimator now uses the same treated units as
  the fit. A treated unit must have a treated cell with an outcome. Before, a
  unit with `W = 1` but no treated outcome changed the weights of each draw,
  also in version 0.3.0. Thus, a draw with the rows of the data gave a
  different estimate.
- The global placebo cross-validation now uses the placebo sets that it
  draws. Before, with `placebo_size` more than 1, the C code read the units of
  the sets in the wrong order. Then a set sometimes had the same unit two
  times. This also occurred in version 0.3.0.
- In the global placebo cross-validation, a placebo unit is now a treated unit
  only if it has an outcome in the held-out periods. The fit on the data uses
  the same rule. Before, a placebo unit without such an outcome changed the
  weights. This also occurred in version 0.3.0.
- The mean of very large effects is now finite. The RMSE of the global placebo
  cross-validation is also finite for very large placebo effects. Before, the
  sums were sometimes infinite.

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
