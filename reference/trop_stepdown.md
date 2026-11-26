# Romano-Wolf stepdown p-values

`trop_stepdown()` adjusts the p-values of many tests for multiple
testing. The tests can be for many outcomes, or for many groups of one
outcome. The function uses the stepdown method of Romano and Wolf
(2005). It controls the probability of one or more false rejections. The
method uses the joint bootstrap distribution of all estimates. Thus, all
estimates must have the same bootstrap draws.

## Usage

``` r
trop_stepdown(x)
```

## Arguments

- x:

  A `trop` object, a `trop_effects` object with a bootstrap, or a list
  of such objects. Use a named list to name the tests. All objects must
  have the same bootstrap rows. To make them, give the same `boot_rows`
  from
  [`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md)
  to each
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) and
  [`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
  call.

## Value

A data frame with one row for each estimate and with the columns `name`,
`estimate`, `std.error`, `t`, `p_raw`, and `p_adj`. The column `p_raw`
is the bootstrap p-value of one test, and `p_adj` is the adjusted
p-value. An estimate with no standard error has no p-value.

## Details

For each estimate, the function calculates the statistic \\t =
\hat\theta / se\\, where the standard error is the standard deviation of
the draws. The null hypothesis is that the effect is 0. The statistics
of the draws are the absolute values of the centered draws divided by
the standard errors. The function sorts the estimates by \\\|t\|\\. In
each step, it compares \\\|t\|\\ with the largest statistic of the draws
over the estimates that remain.

The p-value is (1 + the number of draws with a larger value) / (the
number of draws + 1). An adjusted p-value is never smaller than the
adjusted p-value of an estimate with a larger \\\|t\|\\.

The function uses only the draws that have a value for all estimates.
The result is not exact if the number of draws is small.

## References

Romano, J. P., and Wolf, M. (2005). Stepwise Multiple Testing as
Formalized Data Snooping. *Econometrica*, 73(4), 1237-1282.

## See also

[`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md),
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                     n_post = 3, effect = 1, seed = 1)
W <- tapply(sim$d, list(sim$unit, sim$time), sum)
rows <- trop_boot_rows(W, n_boot = 40, seed = 1)
fit <- trop(sim, "y", "d", "unit", "time",
            lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
            boot_rows = rows, control = trop_control(threads = 1))
trop_stepdown(trop_effects(fit, by = "event_time"))
#>   name    estimate std.error           t      p_raw      p_adj
#> 1    0  1.41712290 0.3612847  3.92245474 0.02439024 0.02439024
#> 2    1  1.97684134 0.4881943  4.04929179 0.02439024 0.02439024
#> 3    2 -0.03762763 0.3831373 -0.09820927 0.87804878 0.87804878
```
