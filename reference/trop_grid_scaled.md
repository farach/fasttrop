# Scale the tuning grid to the data

`trop_grid_scaled()` changes the values of a tuning grid so that they
fit the scale of the outcome. The default grid of
[`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md)
uses fixed values. These values fit an outcome with a standard deviation
of about 1. For an outcome with other units, the best value of a
parameter can be outside of the grid. Then cross-validation chooses the
edge of the grid.

## Usage

``` r
trop_grid_scaled(Y, W, grid = trop_grid(), demean = FALSE, max_units = 500L)
```

## Arguments

- Y:

  A numeric matrix of outcomes with units in rows. `NA` shows a cell
  without an outcome.

- W:

  A matrix of treatment indicators with the values 0 and 1. It must have
  the same dimensions as `Y`.

- grid:

  The tuning grid from
  [`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md)
  to scale.

- demean:

  Logical. If `TRUE`, the function subtracts the mean of the untreated
  outcomes of each unit before it calculates the scales. Refer to the
  argument `demean` of
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

- max_units:

  The maximum number of units for the median gap.

## Value

A list of class `trop_grid` with the attribute `scales`. The attribute
is a list with the elements `unit` (the median gap) and `nn` (the
residual scale).

## Details

The function changes the grid in two steps:

- The unit weights decrease with the root mean square gap between the
  outcome histories of two units. The function divides the values of
  `unit` by the median of this gap. The gap uses the untreated cells
  that both units have. For a panel with more than `max_units` units,
  the median uses `max_units` units that are the same distance apart in
  the order of the rows.

- The nuclear-norm penalty depends on the scale of the noise. The
  function multiplies the finite values of `nn` by a scale. The scale is
  the root mean square of the residuals of a two-way fixed effects fit
  to the untreated cells.

The values of `time` and the value `Inf` do not change. The scales are
in the attribute `scales` of the result.

With `grid = "auto"`,
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) and
[`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
use this function with the default grid, with the same `demean` and the
same `target`.

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                     n_post = 3, effect = 1, seed = 1)
pm <- tapply(sim$y, list(sim$unit, sim$time), sum)
W <- tapply(sim$d, list(sim$unit, sim$time), sum)
g <- trop_grid_scaled(pm * 100, W)
attr(g, "scales")
#> $unit
#> [1] 165.3735
#> 
#> $nn
#> [1] 92.0611
#> 
```
