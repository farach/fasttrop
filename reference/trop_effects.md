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
  seed = NULL
)

# S3 method for class 'trop_effects'
print(x, digits = 4L, ...)

# S3 method for class 'trop_effects'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
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

- x:

  A `trop_effects` object.

- digits:

  The number of significant digits to show.

- ...:

  Not used.

- row.names, optional:

  Not used.

## Value

A data frame of class `trop_effects` with one row for each group. The
columns are the group (the name of the column is the value of `by`),
`estimate`, and `n_cells`. With a bootstrap, the columns `std.error`,
`conf.low`, `conf.high`, `band.low`, and `band.high` are also in the
data frame. The data frame has these attributes: `draws` (a matrix with
one column for each group), `rows` (the rows of the draws), and
`critical` (the critical value of the band).

## Details

The function uses the effects `tau` of the fit. Thus, the target cells
of the fit are the only cells in the groups. A fit that has an `NA`
effect for a cell does not give a value for that cell.

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

The argument `window` limits the cells to cells with an event time in
the window. For example, `window = c(0, 3)` uses the first 4 treated
periods.

## Inference

The bootstrap uses the unit draws of the fit. For each draw, the
function calculates the local estimator again with the tuning parameters
of the fit. Then it averages the effects of each group. The tuning
parameters do not change. The standard error of a group is the standard
deviation of its draws. A group that has no treated cell in a draw has
no value for that draw.

By default, the function uses the draws that the fit stored in
`fit$boot$rows`. Then the average of all treated cells has the same
draws as `fit$boot$draws`. With `n_boot`, the function makes new draws.
With `boot_rows`, you give the draws. A group of one fit and a group of
another fit have a joint distribution only if they use the same draws.

The columns `band.low` and `band.high` give a simultaneous band with the
sup-t method. The band has a critical value for all groups. The critical
value is a quantile of the largest absolute value of the standardized
draws over the groups. The band contains all true group effects with the
probability `conf_level`, if the bootstrap is correct.

The intervals `conf.low` and `conf.high` are for one group each. The
band needs many draws. With only one group, the band columns are `NA`.

## See also

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
#>  event_time estimate n_cells std.error conf.low conf.high band.low band.high
#>           0  1.41712       5    0.3535   0.7243     2.110   0.6428    2.1914
#>           1  1.97684       5    0.4643   1.0669     2.887   0.9599    2.9938
#>           2 -0.03763       5    0.3825  -0.7873     0.712  -0.8755    0.8002
#> Bootstrap: 50 draws. The band is simultaneous (sup-t) at the 95% level.
```
