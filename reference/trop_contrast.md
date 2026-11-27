# Linear combinations of group effects

`trop_contrast()` calculates linear combinations of the estimates of
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
or [`trop()`](https://farach.github.io/fasttrop/reference/trop.md)
objects. Examples are the mean effect of event times 0 to 3, or the
difference between two cohorts. The function uses the bootstrap draws of
the objects. It does not calculate the estimator again.

## Usage

``` r
trop_contrast(x, contrasts, conf_level = NULL, ci = c("normal", "percentile"))
```

## Arguments

- x:

  A `trop_effects` object, a `trop` object, or a named list of such
  objects, all with a bootstrap.

- contrasts:

  A named numeric vector of coefficients, with the names of the
  estimates. A named list of such vectors gives one combination for each
  element. An estimate that is not in a vector has the coefficient 0.

- conf_level:

  The confidence level. The default is the level of the first object.

- ci:

  Use `"normal"` or `"percentile"`. Refer to
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

## Value

A data frame of class `trop_effects` with one row for each combination
and the columns `contrast`, `estimate`, `std.error`, `conf.low`,
`conf.high`, `band.low`, `band.high`, and `n_draws`. The attribute
`coefficients` is the matrix of coefficients.

## Details

For each combination, the estimate is the sum of the coefficients times
the estimates. Each draw gets the same sum of its draws. A draw without
a value for an estimate with a coefficient that is not 0 has no value
for the combination. The standard error is the standard deviation of the
draws.

All objects must have the same bootstrap draws. Draw b of each object
must use the same units. Refer to the section Joint draws of
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md).
The names of the estimates are the names of
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md).
For one `trop_effects` object, the names are the groups, for example
`"0"` and `"1"` for event times. For a named list, the names have the
form `"event: 0"`.

The result is a `trop_effects` object. Thus,
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)
can adjust the p-values of the combinations, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) gives their covariance.

## See also

[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                     n_post = 3, effect = 1, seed = 1)
fit <- trop(sim, "y", "d", "unit", "time",
            lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
            n_boot = 50, seed = 1, control = trop_control(threads = 1))
ev <- trop_effects(fit, by = "event_time")
trop_contrast(ev, list(mean_0_2 = c("0" = 1, "1" = 1, "2" = 1) / 3,
                       last_minus_first = c("2" = 1, "0" = -1)))
#> TROP effects by contrast
#>   Each row is a linear combination of estimates.
#>          contrast estimate std.error conf.low conf.high band.low band.high
#>          mean_0_2    1.119    0.2020   0.7228    1.5147   0.7323    1.5052
#>  last_minus_first   -1.455    0.5742  -2.5801   -0.3294  -2.5531   -0.3564
#>  n_draws
#>       50
#>       50
#> Bootstrap: 50 draws. The band is simultaneous (sup-t) at the 95% level. It uses the 50 draws with a value for all rows.
```
