# Choose a tuning grid

Suppose you change an outcome from dollars to cents. The units and
periods stay the same, but each gap between two outcomes becomes 100
times as large. The default grid still tests the same values. Use
`grid = "auto"` to adjust the grid to the scale of the untreated
outcomes.

## Fit a practice panel

This example uses a panel with one treatment start date.
Cross-validation compares settings on a sample of untreated cells. The
example omits the bootstrap because the question here is which settings
to use.

``` r

library(fasttrop)
panel <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                       n_post = 3, seed = 1)
panel$cents <- 100 * panel$y
fit <- trop(panel, "cents", "d", "unit", "time", grid = "auto",
            cv_cells = 30, se = "none", seed = 1)
fit$lambda
#>         time         unit           nn 
#>  0.500000000  0.003023458 92.061103545
fit$cv$boundary
#> NULL
```

`fit$lambda` gives the values that the search chose. If
`fit$cv$boundary` names a parameter, the search chose a value at one end
of its range. This is not proof that the fit failed. It means that the
search did not compare values beyond that end.

## Change the range

If the warning names the largest value of `unit`, add a larger value to
that part of the grid. Keep the other values for comparison. The code
below starts with the scaled grid and adds one candidate. Use the same
seed and sample size to compare settings on the same cells.

``` r

Y <- tapply(panel$cents, list(panel$unit, panel$time), sum)
W <- tapply(panel$d, list(panel$unit, panel$time), sum)
g <- trop_grid_scaled(Y, W)
wider <- trop_grid(time = g$time, unit = c(g$unit, 2 * max(g$unit)),
                   nn = g$nn)
refit <- trop(panel, "cents", "d", "unit", "time", grid = wider,
              cv_cells = 30, se = "none", seed = 1)
rbind(first = fit$lambda, wider = refit$lambda)
#>       time        unit      nn
#> first  0.5 0.003023458 92.0611
#> wider  0.5 0.003023458 92.0611
```

For a warning at the smallest positive value, try a smaller positive
value instead. Choose the range from the prediction errors on untreated
cells, not from the size of the treatment effect.

## What the scale changes

[`trop_grid_scaled()`](https://farach.github.io/fasttrop/reference/trop_grid_scaled.md)
divides the `unit` values by a typical gap between untreated histories.
It multiplies finite `nn` values by the residual scale after it removes
unit and period effects. It does not change `time` or `Inf`. To inspect
the two scales, use `attr(g, "scales")`.

The scaled grid does not remove differences in level between units. If
you want distances to reflect changes rather than levels, set
`demean = TRUE`. The fit then subtracts each unit’s mean untreated
outcome. When you supply a scaled grid yourself, use the same `demean`
setting for the grid and the fit.
