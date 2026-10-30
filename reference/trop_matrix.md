# TROP estimator for outcome and treatment matrices

`trop_matrix()` is the matrix version of
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md). The
rows are units and the columns are periods.

## Usage

``` r
trop_matrix(
  Y,
  W,
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

- Y:

  A numeric matrix of outcomes with one row for each unit and one column
  for each period. `NA` shows a cell without an outcome.

- W:

  A matrix of treatment indicators with the values 0 and 1. It must have
  the same dimensions as `Y`.

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

An object of class `trop`. Refer to
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(15 * 10), 15, 10)
W <- matrix(0L, 15, 10)
W[1:3, 8:10] <- 1L
fit <- trop_matrix(Y, W, lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
                   control = trop_control(threads = 1))
coef(fit)
#>      ATT 
#> 0.110928 
```
