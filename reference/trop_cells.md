# Count the cells of each group before the bootstrap

`trop_cells()` shows which treated cells are in each group of
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
and which cells are not. It does not show the effects. Thus, you can
examine the estimand before you see the results or calculate the
bootstrap.

## Usage

``` r
trop_cells(
  fit,
  by = c("event_time", "cohort", "time", "unit", "overall"),
  window = NULL,
  anticipation = 0L,
  weights = c("cell", "unit"),
  unit_weights = NULL
)

# S3 method for class 'trop_cells'
print(x, digits = 3L, ...)
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

- x:

  A `trop_cells` object.

- digits:

  The number of significant digits to show.

- ...:

  Not used.

## Value

A data frame of class `trop_cells` with one row for each group and the
columns `n_units`, `n_cells`, `n_excluded`, `n_skipped`, and `share`.
The column `n_excluded` counts the treated cells with an outcome that
are outside of the target. The column `n_skipped` counts the cells in
the target that have no effect, because the untreated cells do not
identify them or their fit failed. The attribute `totals` has the counts
for all treated cells.

## Details

The function uses the groups, the window, the anticipation, and the
weights of
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md).
For each group, it counts the units and the cells in the estimate. It
also counts the treated cells outside of the target of the fit, and the
cells that the untreated cells do not identify. A group with only cells
outside of the estimate is also in the result.

The column `share` is the weight of each group in the overall estimate:
the mean of all cells in the window with the same weights. With equal
weights for each cell, the share of a group is its number of cells
divided by the number of cells in the estimate.

If no treated cell is in a group, the result has no rows. The attribute
`totals` and the print method then show where the treated cells are.

## See also

[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 6,
                     n_post = 4, design = "staggered", seed = 1)
sim$in_target <- sim$time <= 11
fit <- suppressWarnings(
  trop(sim, "y", "d", "unit", "time", lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
       target = "in_target", se = "none", control = trop_control(threads = 1))
)
trop_cells(fit, by = "event_time")
#> TROP cells by event time
#>   Weights: each cell has the same weight.
#>  event_time n_units n_cells n_excluded n_skipped share
#>           0       4       4          2         0 0.667
#>           1       2       2          2         0 0.333
#>           2       0       0          2         0 0.000
#> Cells in the estimate: 6. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 6 outside of the target.
```
