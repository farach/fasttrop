# Triply Robust Panel (TROP) estimator

`trop()` calculates the average treatment effect on the treated (ATT)
with the TROP estimator of Athey, Imbens, Qu, and Viviano (2026). A cell
is one unit in one period. For each treated cell, a model calculates the
outcome without the treatment. The model has two-way fixed effects and a
low-rank term with a nuclear-norm penalty, and it uses only untreated
cells. Unit weights decrease with the distance between the outcome
history of a unit and the treated unit. Time weights decrease with the
distance in periods.

## Usage

``` r
trop(
  data,
  outcome,
  treatment,
  unit,
  time,
  method = c("local", "global"),
  lambda = NULL,
  grid = trop_grid(),
  cv_cells = NULL,
  n_placebo = 200L,
  placebo_size = 1L,
  max_cycles = 10L,
  se = c("bootstrap", "none"),
  n_boot = 200L,
  resample = c("stratified", "controls"),
  conf_level = 0.95,
  ci = c("normal", "percentile"),
  control = trop_control(),
  seed = NULL
)
```

## Arguments

- data:

  A data frame in long format with one row for each unit and period.

- outcome, treatment, unit, time:

  Column names as strings. The treatment must be 0 or 1. Outcomes can be
  `NA`. The fits do not use these cells.

- method:

  Use `"local"` or `"global"`. Refer to the section Methods.

- lambda:

  The tuning parameters as a named numeric vector or a named list with
  the names `time`, `unit`, and `nn`. If `lambda` is `NULL`,
  cross-validation selects the tuning parameters from `grid`.

- grid:

  The tuning grids from
  [`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md).

- cv_cells:

  For the local method, the number of untreated cells in a random sample
  for the leave-one-out criterion. If `cv_cells` is `NULL`, the
  criterion uses all untreated cells, as in Eq. 5 of the paper.

- n_placebo, placebo_size:

  For the global method, the number of placebo sets and the number of
  untreated units in each set. The defaults are the values of the Python
  package from the authors: 200 sets of 1 unit.

- max_cycles:

  The maximum number of cycles of the coordinate search.

- se:

  Use `"bootstrap"` for bootstrap standard errors or `"none"` for no
  standard errors.

- n_boot:

  The number of bootstrap draws.

- resample:

  Use `"stratified"` or `"controls"`. Refer to the section Inference.

- conf_level:

  The confidence level of the interval.

- ci:

  Use `"normal"` or `"percentile"`. The normal interval is the estimate
  plus or minus a normal quantile times the bootstrap standard error.
  The percentile interval uses quantiles of the bootstrap draws. Refer
  to the section Inference.

- control:

  The solver options from
  [`trop_control()`](https://farach.github.io/fasttrop/reference/trop_control.md).

- seed:

  An optional integer seed for the sampled cells, the placebo sets, and
  the bootstrap. The function does not change the state of the random
  number generator in the R session.

## Value

An object of class `trop`. It is a list with these elements and with
other elements:

- `estimate`: the ATT.

- `std.error`, `conf.low`, and `conf.high`: the bootstrap standard error
  and the confidence interval.

- `tau`: a matrix with the effect for each treated cell.

- `lambda`: the tuning parameters.

- `cv`: the results of the cross-validation, or `NULL` for fixed tuning
  parameters.

- `boot`: the bootstrap draws.

- `timing`: the time in the C code and the total time, in seconds.

## Methods

- The local method (`method = "local"`) is the estimator of Eq. 12 and
  Algorithm 2 of the paper. Each treated cell has its own unit weights,
  time weights, and fit. The leave-one-out criterion of Eq. 4 and Eq. 5
  selects the tuning parameters. The coordinate search follows footnote
  2 of the paper. You can use this method for block, staggered, and
  general designs.

- The global method (`method = "global"`) is the estimator of the Python
  package from the authors. It calculates one fit for a block design.
  The time weights decrease with the distance from the center of the
  treated block. The unit weights decrease with the distance from the
  average treated unit before the treatment. Placebo cross-validation on
  untreated units selects the tuning parameters. The coordinate search
  starts at the means of the grids.

## Designs

In a block design, all treated units start the treatment in the same
period and stay treated to the end of the panel. In a staggered design,
the start periods are different, but the treatment does not stop. In a
general design, the treatment stops for one or more units.

## Inference

The unit bootstrap of Algorithm 3 of the paper gives the standard error.
The bootstrap keeps the tuning parameters fixed. With
`resample = "stratified"`, the bootstrap draws treated and untreated
units separately. With `resample = "controls"`, the treated units do not
change. The standard error is the standard deviation of the draws.

The default interval is the normal interval with this standard error, as
in the bootstrap variance of the paper. With `ci = "percentile"`, the
interval comes from quantiles of the draws, as in the Python package
from the authors. A percentile interval is less accurate with a small
number of draws. The paper also says that the bootstrap is correct only
with many treated units. Thus, the function gives a warning message for
a bootstrap with fewer than 5 treated units.

## Limits of the triple robustness result

Theorem 5.1 of the paper gives a bound for the bias. The bound is the
product of three terms: the unit imbalance, the time imbalance, and the
bias of the regression adjustment. The proof uses fixed weights, and it
is only for a block design. For a staggered or general design, the
function calculates the estimator, but the bound has no proof. Thus, the
function gives a warning message for these designs.

## References

Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## Examples

``` r
sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
                     n_post = 3, effect = 1, seed = 1)
fit <- trop(sim, "y", "d", "unit", "time",
            lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
            n_boot = 50, seed = 1, control = trop_control(threads = 1))
fit
#> Triply Robust Panel (TROP) estimate
#>   Method: local | Design: block | Units: 20 | Periods: 12 | Treated cells: 15
#>   ATT: 1.119  (SE 0.202, 95% CI 0.7228 to 1.515)
#>   lambda: time = 0.2, unit = 0.5, nn = 0.5
```
