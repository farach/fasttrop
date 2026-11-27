# Effects by event time, cohort, period, or unit

`trop_effects()` averages the effects of the treated cells of a
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) fit in
groups. A group is, for example, all treated cells with the same number
of periods since the start of the treatment. The function can also give
bootstrap standard errors, intervals, and simultaneous bands for all
groups.

## Usage

``` r
trop_effects(
  fit,
  by = c("event_time", "cohort", "time", "unit", "overall"),
  window = NULL,
  n_boot = NULL,
  boot_rows = NULL,
  strata = NULL,
  resample = NULL,
  conf_level = fit$conf.level,
  ci = c("normal", "percentile"),
  seed = NULL,
  anticipation = 0L,
  weights = c("cell", "unit"),
  unit_weights = NULL,
  checkpoint = NULL,
  progress = FALSE
)

# S3 method for class 'trop_effects'
print(x, digits = 4L, ...)

# S3 method for class 'trop_effects'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'trop_effects'
coef(object, ...)

# S3 method for class 'trop_effects'
vcov(object, ...)
```

## Arguments

- fit:

  A `trop` object from
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) or
  [`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
  with `method = "local"`.

- by:

  The groups. Refer to the section Groups.

- window:

  An optional vector with 2 event times. The function uses cells with an
  event time from `window[1]` to `window[2]`.

- n_boot:

  The number of new bootstrap draws. With `NULL`, the function uses the
  draws of the fit, if the fit has draws. With 0, the function does not
  calculate a bootstrap.

- boot_rows:

  A matrix of draws from
  [`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md).
  This argument replaces `n_boot`.

- strata, resample:

  The settings for new draws. Refer to
  [`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md).
  With `resample = NULL`, the function uses the setting of the fit.

- conf_level:

  The confidence level. The default is the level of the fit.

- ci:

  Use `"normal"` or `"percentile"`. Refer to
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

- seed:

  An optional integer seed for new draws. The function does not change
  the state of the random number generator in the R session.

- anticipation:

  The number of treated periods at the start of each run that come
  before the true start. Refer to the section Anticipation. The default
  0 uses the first treated period as the start.

- weights:

  Use `"cell"` or `"unit"`. Refer to the section Weights.

- unit_weights:

  An optional numeric vector with one weight of 0 or more for each unit.
  A named vector must have a value for each unit name of the fit. A
  vector without names must be in the order of the rows of the panel.

- checkpoint:

  An optional path of a file for the bootstrap draws. Refer to the
  section Long runs.

- progress:

  Logical. If `TRUE`, the function shows the progress of the bootstrap.

- x:

  A `trop_effects` object.

- digits:

  The number of significant digits to show.

- ...:

  Not used.

- row.names, optional:

  Not used.

- object:

  A `trop_effects` object.

## Value

A data frame of class `trop_effects` with one row for each group. The
columns are the group (the name of the column is the value of `by`),
`estimate`, `n_cells`, and `n_units`. With a bootstrap, the columns
`std.error`, `conf.low`, `conf.high`, `band.low`, `band.high`, and
`n_draws` are also in the data frame. The data frame has these
attributes: `draws` (a matrix with one column for each group), `rows`
(the rows of the draws), `critical` (the critical value of the band),
`n_complete` (the number of draws with a value for all groups),
`diagnostics`, and `units` (the names of the units of the fit). When
`by` is a list, the result is a named list of such data frames.

[`coef()`](https://rdrr.io/r/stats/coef.html) gives the estimates with
the group names. [`vcov()`](https://rdrr.io/r/stats/vcov.html) gives the
covariance matrix of the draws, from the draws with a value for all
groups.

## Details

The function uses the effects `tau` of the fit. Thus, the target cells
of the fit are the only cells in the groups. A cell that has no effect
in the fit, for example a cell outside of the target, is not in a group.
A group with no such cell is not in the result. Use
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
to count the cells of each group before you calculate the bootstrap.

## Groups

The argument `by` sets the groups:

- `"event_time"`: the number of periods since the start of the current
  treatment. The value 0 is the first treated period. If a unit stops
  the treatment and starts it again, the count starts again.

- `"cohort"`: the period in which the current treatment started.

- `"time"`: the period.

- `"unit"`: the unit.

- `"overall"`: one group with all cells. The result is the ATT.

- A matrix with the dimensions of the panel. The values are the group of
  each cell. A cell with `NA` is not in a group. Only treated cells with
  an effect are in groups.

- A named list of the values above, for example
  `list(event = "event_time", cohort = "cohort")`. The function then
  gives a named list with one result for each element. All results use
  the same bootstrap draws, and the function calculates the estimator
  only one time for each draw. A character element without a name uses
  its value as the name.

The argument `window` limits the cells to cells with an event time in
the window. For example, `window = c(0, 3)` uses the first 4 treated
periods.

## Anticipation

Sometimes the treatment has an effect before it starts, for example when
units know the start date. A usual procedure is to set `W = 1` in the
periods before the start, so that these cells are not controls. Then,
the first `W = 1` cell is not the true start.

With `anticipation = k`, the function counts the first `k` treated
periods of each run as the periods before the start. Thus, these periods
get the event times `-k` to `-1`, and the event time 0 is the true
start. The cohort is also the period of the true start. A cohort that
starts after the last period of the panel is not in a group.

To keep the anticipation cells out of the estimate, also remove them
from the `target` of the fit. Use
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)
to make `W` and `target` from the start dates.

The function finds the start of each run from `W`. If a run starts in
the first period of the panel, its anticipation periods can be before
the panel. Then the event times of the run can be wrong, and the
function gives a warning message. If the fit has event times from the
argument `event_time` of
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md), the
function uses them. Then `anticipation` must be 0. Refer to the section
Event times of
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

## Weights

By default, each cell in a group has the same weight. Thus, a unit with
more cells in a group has more weight in the group. With
`weights = "unit"`, each unit has the same total weight in each group.
The function divides the weight of a unit by its number of cells in the
group.

The weights are relative. In each group, the function divides the
weights by the largest weight of the units of that group. Thus, the
estimate of a group does not depend on the weights of other units. In a
group, the largest positive weight must be less than 1e200 times the
smallest positive weight.

With `unit_weights`, each unit has its own weight, for example its
number of workers before the treatment. Use weights that the treatment
cannot change. With `weights = "cell"`, each cell gets the weight of its
unit. With `weights = "unit"`, the weight of a unit is its total weight
in each group. A unit with a weight of 0 is not in a group.

The bootstrap uses the same weights. A unit that a draw selects two
times counts as two units in the draw. The ATT of the fit always gives
the same weight to each cell.

## Inference

The bootstrap uses the unit draws of the fit. For each draw, the
function calculates the local estimator again with the tuning parameters
of the fit. Then it averages the effects of each group. The tuning
parameters do not change. The standard error of a group is the standard
deviation of its draws.

A group that has no treated cell in a draw has no value for that draw.
The column `n_draws` gives the number of draws with a value for each
group.

The draws run in parallel with the `threads` of the control of the fit.

By default, the function uses the draws that the fit stored in
`fit$boot$rows`. Then the average of all treated cells has the same
draws as `fit$boot$draws`. With `n_boot`, the function makes new draws.
With `boot_rows`, you give the draws. A group of one fit and a group of
another fit have a joint distribution only if they use the same draws.

Each call calculates the estimator again for each draw. To get many sets
of groups, give them in one list to `by`. To calculate the draws only
one time, make the fit with `se = "none"`. Then include `"overall"` in
the list, and give the draws with `boot_rows`.

The columns `band.low` and `band.high` give a simultaneous band with the
sup-t method. The band has a critical value for all groups. The critical
value is a quantile of the largest absolute value of the standardized
draws over the groups. The band contains all true group effects with the
probability `conf_level`, if the bootstrap is correct.

The band uses only the draws with a value for all groups. The centers
and the standard errors of the band also come from these draws. Thus, if
some draws have no value for a group, the band and the interval of that
group use different draws. If fewer than 2 draws have a value for all
groups, the band is not available.

The intervals `conf.low` and `conf.high` are for one group each. The
band needs many draws. With only one group, the band columns are `NA`.

## Diagnostics

The attribute `diagnostics` is a data frame with one row for each draw.
Its columns count these cells and fits of the draw:

- `n_used`: the treated cells with an effect.

- `n_skipped`: the cells that the untreated cells do not identify.

- `n_nonconverged`: the fits that did not converge.

- `n_failed`: the fits that failed. A fit fails if the singular value
  decomposition fails or if the result is not finite.

- `n_dropped`: the cells that have an effect in the fit, but that the
  untreated cells of the draw do not identify. The count does not
  include the fits that failed.

The function gives a warning message when a draw has such fits or cells.
It also gives a warning message for a group with fewer than 2 draws with
a value. The print method shows these counts.

## Long runs

With `progress = TRUE`, the function shows the number of draws that are
done and an estimate of the time that is left. With `checkpoint`, the
function saves the draws to a file after each block of draws. If the
calculation stops, call the function again with the same arguments and
the same file. The function then calculates only the draws that are not
in the file. The results are equal to the results of one run without a
stop.

The file has a fingerprint of the calculation. The fingerprint uses the
panel, the tuning parameters, the solver settings, the draw rows, the
groups, the weights, and the version of fasttrop. If the fingerprint of
the file is different, the function stops with an error and does not
change the file. The number of threads is not in the fingerprint,
because it does not change the results.

The file is an RDS file with the draws, but not the panel. Keep it in a
folder on your computer. The function makes sure that it can write to
the folder before it starts the bootstrap.

## See also

[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md),
[`trop_contrast()`](https://farach.github.io/fasttrop/reference/trop_contrast.md),
[`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md),
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                     n_post = 3, effect = 1, seed = 1)
fit <- trop(sim, "y", "d", "unit", "time",
            lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
            n_boot = 50, seed = 1, control = trop_control(threads = 1))
trop_effects(fit, by = "event_time")
#> TROP effects by event time
#>   Weights: each cell has the same weight.
#>  event_time estimate n_cells n_units std.error conf.low conf.high band.low
#>           0  1.41712       5       5    0.3535   0.7243     2.110   0.6428
#>           1  1.97684       5       5    0.4643   1.0669     2.887   0.9599
#>           2 -0.03763       5       5    0.3825  -0.7873     0.712  -0.8755
#>  band.high n_draws
#>     2.1914      50
#>     2.9938      50
#>     0.8002      50
#> Bootstrap: 50 draws. The band is simultaneous (sup-t) at the 95% level. It uses the 50 draws with a value for all rows.

# Two sets of groups with one bootstrap sweep.
both <- trop_effects(fit, by = list(event = "event_time", time = "time"))
trop_stepdown(both)
#> Romano-Wolf stepdown p-values
#>      name estimate std.error        t   p_raw   p_adj
#>  event: 0  1.41712    0.3535  4.00906 0.01961 0.01961
#>  event: 1  1.97684    0.4643  4.25805 0.01961 0.01961
#>  event: 2 -0.03763    0.3825 -0.09838 0.90196 0.90196
#>  time: 10  1.41712    0.3535  4.00906 0.01961 0.01961
#>  time: 11  1.97684    0.4643  4.25805 0.01961 0.01961
#>  time: 12 -0.03763    0.3825 -0.09838 0.90196 0.90196
#> Draws: 50 of 50 have a value for all estimates. The p-values use these draws.
```
