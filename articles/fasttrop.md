# Get started with fasttrop

This vignette shows a first analysis with fasttrop. You make a practice
panel, you calculate the TROP estimate, and you read the result. The
article “TROP in plain words” on the [package
website](https://farach.github.io/fasttrop/) explains the ideas of the
method.

## The data

fasttrop uses panel data in long format. Each row is one unit in one
period. The data must have four columns:

- the outcome
- the treatment indicator, with the values 0 and 1
- the unit identifier
- the period identifier

The function
[`trop_simulate()`](https://farach.github.io/fasttrop/reference/trop_simulate.md)
makes a practice panel with a known effect. In this panel, 5 of 20 units
get the treatment in the last 3 of 12 periods. The true effect is 1 for
each treated cell.

``` r

library(fasttrop)
panel <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5, n_post = 3,
                       effect = 1, seed = 1)
head(panel)
#>   unit time          y d         y0
#> 1    1    1 -0.5496887 0 -0.5496887
#> 2    2    1 -0.5309372 0 -0.5309372
#> 3    3    1 -1.1621093 0 -1.1621093
#> 4    4    1  2.5570765 0  2.5570765
#> 5    5    1 -0.1638954 0 -0.1638954
#> 6    6    1 -1.1801013 0 -1.1801013
```

The column `y` is the outcome, and the column `d` is the treatment
indicator. The column `y0` is the outcome without the treatment. Real
data does not have this column. The practice panel has it, so you can
compare the estimate with the truth.

## A first estimate

To calculate the estimate, give
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) the data
and the names of the four columns. In this first fit, you set the tuning
parameters. This makes the fit fast.

``` r

fit <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
            lambda = c(time = 0.2, unit = 0.5, nn = 0.5), n_boot = 100, seed = 1)
fit
#> Triply Robust Panel (TROP) estimate
#>   Method: local | Design: block | Units: 20 | Periods: 12 | Treated cells: 15
#>   ATT: 1.119  (SE 0.2199, 95% CI 0.6879 to 1.55)
#>   lambda: time = 0.2, unit = 0.5, nn = 0.5
```

After the title, the output has three lines:

1.  The method, the design of the panel, the size of the panel, and the
    number of treated cells.
2.  The average treatment effect on the treated (ATT), with the
    bootstrap standard error (SE) and the 95% confidence interval.
3.  The tuning parameters.

The true ATT of this panel is 1. The interval contains this value.

## Parts of the result

Use these functions to get parts of the result:

``` r

coef(fit)
#>      ATT 
#> 1.118779
confint(fit)
#>         2.5 %   97.5 %
#> ATT 0.6878718 1.549686
head(as.data.frame(fit))
#>   unit time        tau
#> 1    7   10  2.0694817
#> 2    7   11  3.3467710
#> 3    7   12  0.1396825
#> 4   11   10  0.7607983
#> 5   11   11  1.9304028
#> 6   11   12 -0.3167042
```

The function
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
effect for each treated cell. The ATT is the mean of these effects.

## Let the data select the tuning parameters

TROP has three tuning parameters:

- `time` sets how fast the time weights decrease with the distance from
  the treated period. With 0, all periods get the same weight.
- `unit` sets how fast the unit weights decrease when the outcome
  history of a unit is different from the treated unit. With 0, all
  units get the same weight.
- `nn` is the penalty on the low-rank term. A larger value gives a
  simpler low-rank term. With `Inf`, the model has no low-rank term.

If you do not set `lambda`, cross-validation selects the tuning
parameters from the grid of
[`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md).
To make this example fast, the cross-validation uses a sample of 60
untreated cells.

``` r

fit_cv <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
               cv_cells = 60, n_boot = 100, seed = 1)
summary(fit_cv)
#> Triply Robust Panel (TROP) estimate
#>   Method: local | Design: block | Units: 20 | Periods: 12 | Treated cells: 15
#>   ATT: 1.142  (SE 0.2267, 95% CI 0.6978 to 1.587)
#>   lambda: time = 0.1, unit = 0.5, nn = 1
#>   Tuning: loocv, score 18.79, 34 grid evaluations
#>   Bootstrap: 100 draws (stratified), 0 failed, normal interval
#>   Treated units: 5 | Skipped cells: 0 | Non-converged fits: 0
#>   Time: 1.130 s in total, 1.126 s in the C code (4 threads)
```

The function [`summary()`](https://rdrr.io/r/base/summary.html) also
gives the result of the tuning, the bootstrap, and the time of the
calculation. With `cv_cells = NULL`, the default, the cross-validation
uses all untreated cells, as in the paper. This takes more time, but it
uses all the data.

## The local method and the global method

The default method is the local method of the paper. It calculates a
separate fit for each treated cell. The global method is the method of
the Python package from the authors. It calculates one fit for all
treated cells. You can use the global method only for a block design. In
a block design, all treated units start the treatment in the same period
and stay treated.

``` r

fit_global <- trop(panel, outcome = "y", treatment = "d", unit = "unit",
                   time = "time", method = "global",
                   lambda = c(time = 0.2, unit = 0.5, nn = 0.5), n_boot = 100,
                   seed = 1)
coef(fit_global)
#>      ATT 
#> 1.156408
```

## Standard errors and intervals

A unit bootstrap gives the standard error. Each draw takes a sample of
units with replacement and calculates the estimate again. fasttrop keeps
the tuning parameters fixed in each draw. The argument `n_boot` sets the
number of draws. The default is 200 draws.

The default interval is a normal interval. To get a percentile interval,
set `ci = "percentile"`. The bootstrap is correct only with many treated
units. For a panel with fewer than five treated units, fasttrop gives a
warning message.

## Threads and seeds

fasttrop runs its model fits in C, and it can use more than one thread.
The function
[`trop_threads()`](https://farach.github.io/fasttrop/reference/trop_threads.md)
gives the default number of threads. To set the number of threads, use
[`trop_control()`](https://farach.github.io/fasttrop/reference/trop_control.md):

``` r

fit <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
            control = trop_control(threads = 4))
```

The number of threads does not change the results. With the same `seed`,
a fit gives the same result each time. The `seed` argument does not
change the random number generator of your R session.

## Next steps

- The article “How to read your results” explains each part of the
  output and each warning message.
- The article “What TROP can and cannot do” gives the limits of the
  method.
- The help page
  [`?trop`](https://farach.github.io/fasttrop/reference/trop.md) gives
  all arguments.
