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
  seed = NULL,
  covariates = NULL,
  covariate_method = c("joint", "residualize"),
  target = NULL,
  strata = NULL,
  boot_rows = NULL,
  demean = FALSE
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
  [`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md),
  or `"auto"`. With `"auto"`, the function scales the default grid to
  the data. Refer to the section Tuning grid.

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

- covariates:

  Optional names of numeric covariate columns. Refer to the section
  Covariates.

- covariate_method:

  Use `"joint"` or `"residualize"`. Refer to the section Covariates.

- target:

  An optional name of a column with `TRUE` and `FALSE` (or 1 and 0) for
  each row of `data`. A treated cell with `TRUE` is in the estimate.
  Refer to the section Target cells.
  [`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
  takes a logical or 0/1 matrix with the dimensions of `W`.

- strata:

  An optional name of a column with the stratum of each unit. The value
  must be the same for all rows of a unit.
  [`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
  takes a vector with one value for each unit. Refer to the section
  Inference.

- boot_rows:

  An optional matrix of bootstrap draws from
  [`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md).
  If you use it, `n_boot` and `resample` do not apply, and you cannot
  use `strata`.

- demean:

  Logical. If `TRUE`, the function subtracts the mean of the untreated
  outcomes of each unit before it fits the model. Refer to the section
  Tuning grid.

## Value

An object of class `trop`. It is a list with these elements and with
other elements:

- `estimate`: the ATT.

- `std.error`, `conf.low`, and `conf.high`: the bootstrap standard error
  and the confidence interval.

- `tau`: a matrix with the effect for each treated cell.

- `lambda`: the tuning parameters.

- `cv`: the results of the cross-validation, or `NULL` for fixed tuning
  parameters. The element `cv$boundary` names the tuning values that are
  at the edge of the grid.

- `boot`: the bootstrap draws in `draws`, the number of failed draws in
  `n_failed`, the kind of draws in `resample`, and the draw rows in
  `rows`.

- `n_excluded_cells`: the number of treated cells that are not in the
  target.

- `covariates`: `NULL` without covariates. Otherwise, a list with the
  coefficients `coef`, the `method`, and the names of the covariates
  that the function `used` and `dropped`.

- `panel`: the matrices of the fit.
  [`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
  uses them. They make the object larger.

- `package`: the version of fasttrop, and its commit if it is known.

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

## TWFE imputation

With `lambda = c(time = 0, unit = 0, nn = Inf)`, both methods give TWFE
imputation. This method fits unit effects and period effects to the
untreated cells only. Then it subtracts the fitted value from the
outcome of each treated cell.

A conventional TWFE regression is different. It uses all cells, and it
has one coefficient for the treatment. In a block design with no `NA`
outcomes, the two results are the same. In a staggered design, the
results can be different if the effects are different between treated
cells. The regression can give negative weights to some treated cells.

## Inference

The unit bootstrap of Algorithm 3 of the paper gives the standard error.
The bootstrap keeps the tuning parameters fixed. With
`resample = "stratified"`, the bootstrap draws treated and untreated
units separately. With `resample = "controls"`, the treated units do not
change. The standard error is the standard deviation of the draws.

With `strata`, the bootstrap draws the units of each stratum separately.
Each draw keeps the number of treated units and of untreated units in
each stratum. Use this when the units are in groups, for example a
rollout in waves.

With `boot_rows`, you give the draws. Use
[`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md)
to make them. If you give the same draws to two or more fits, the fits
have a joint bootstrap distribution. Then
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)
can adjust the p-values of the fits. The draws and the panel are in the
result, so
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
can also calculate the effects of groups of cells.

The default interval is the normal interval with this standard error, as
in the bootstrap variance of the paper. With `ci = "percentile"`, the
interval comes from quantiles of the draws, as in the Python package
from the authors. A percentile interval is less accurate with a small
number of draws. The paper also says that the bootstrap is correct only
with many treated units. Thus, the function gives a warning message for
a bootstrap with fewer than 5 treated units.

## Target cells

By default, the estimate is the mean of the effects of all treated
cells. With `target`, the estimate uses only some of the treated cells.
For example, use `target` for the first 3 periods after the start of the
treatment. The function does not use the outcome of a treated cell
outside of the target, because the treatment can change it. The outcome
is not in the fit of any other cell, and the cell is not in the
estimate.

The result has the same value as a fit with `NA` in the outcome of these
cells. The element `n_excluded_cells` of the result has the number of
these cells.

The target is only for the local method, because the global method fits
all treated cells together. Use
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
to calculate the effects of many groups of cells in one fit.

## Tuning grid

The default grid of
[`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md)
fits an outcome with a standard deviation of about 1. If the outcome has
other units, a better value can be outside of the grid. Then
cross-validation chooses a value at the edge of the grid. The function
gives a warning message when this happens for a finite value of `time`,
`unit`, or `nn`. The result has the edge values in `cv$boundary`.

With `grid = "auto"`, the function scales the default grid to the data
with
[`trop_grid_scaled()`](https://farach.github.io/fasttrop/reference/trop_grid_scaled.md).
With `demean = TRUE`, the function subtracts the mean of the untreated
outcomes of each unit before it calculates the distances between units.
The unit effect absorbs a constant of a unit. Thus, with `unit = 0` the
effects do not change. With a larger `unit` weight, the weights of the
units are different. They compare the shape of the outcome histories and
not the level.

## Limits of the triple robustness result

Theorem 5.1 of the paper gives a bound for the bias. The bound is the
product of three terms: the unit imbalance, the time imbalance, and the
bias of the regression adjustment. The proof uses fixed weights, and it
is only for a block design. For a staggered or general design, the
function calculates the estimator, but the bound has no proof. Thus, the
function gives a warning message for these designs.

## Covariates

Section 6.2 of the paper adds a covariate term to the model of Eq. 13.
The model has the fixed effects, the low-rank term, and a linear term
\\X\_{it}'\phi\\. Use `covariates` to give the names of numeric columns.
A row with an outcome must have finite values for all covariates. For a
factor, make indicator columns with
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

The paper gives the model, but it does not give one procedure to
calculate \\\phi\\. The function has two procedures:

- With `covariate_method = "joint"`, each weighted fit calculates
  \\\phi\\ together with the fixed effects and the low-rank term. This
  is the minimization of Eq. 13. The local method calculates one
  \\\phi\\ for each treated cell. The unit weights and the time weights
  use the outcomes, not the outcomes minus \\X'\phi\\. The leave-one-out
  fits and the placebo fits also calculate \\\phi\\.

- With `covariate_method = "residualize"`, the function first calculates
  \\\phi\\ with an unweighted TWFE regression on all untreated cells
  with an outcome. Then it uses \\Y^\* = Y - \tilde{X}\phi\\ for the
  tuning, the weights, and the estimate. In \\\tilde{X}\\, the function
  removes the unit effects and the period effects of the untreated cells
  from \\X\\. The bootstrap calculates \\\phi\\ again for each draw.

Section 6.2 of the paper tunes the estimator after it removes the
covariates from the outcomes. The procedure `"residualize"` is one
interpretation of that text. Usually, it does not give the minimum of
Eq. 13. Its first step uses all untreated outcomes. Thus, the cells that
the cross-validation holds out have an effect on \\\phi\\.

The global method with covariates is an extension of fasttrop. The
Python package from the authors does not have covariates.

The fixed effects absorb a term \\a_i + b_t\\ in a covariate. Thus, the
results do not change when you add a unit constant or a period constant
to a covariate.

The function does not use a covariate that the unit effects, the period
effects, or other covariates explain completely on the untreated cells.
For example, a covariate that does not change in time has no effect on
the estimate. The function gives a warning message for such a covariate.
In a weighted fit with few untreated cells, a covariate can also have no
variation. Then, that fit sets its coefficient to 0.

The element `covariates$coef` of the result contains \\\phi\\ in the
units of the covariates. For the local method with `"joint"`, it is the
mean of the values for the treated cells. This value is only a summary.

The appendix of the paper gives a bias bound with covariates. This bound
adds a term: the imbalance of the covariates times the error in
\\\phi\\. Thus, the weights can decrease the bias, but they do not
remove the bias from a wrong \\\phi\\. A covariate that the treatment
changes can also make a bias. Use only covariates that the treatment
does not change.

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
