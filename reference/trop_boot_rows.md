# Make bootstrap rows

`trop_boot_rows()` makes the matrix of unit indices for the unit
bootstrap of Algorithm 3. Each row of the matrix is one bootstrap draw.
You can give the matrix to the argument `boot_rows` of
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) and
[`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md).
When two or more fits use the same matrix, the draws have the same
units. Then the fits have a joint bootstrap distribution. Use this for a
test of many outcomes, for example with
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md).

## Usage

``` r
trop_boot_rows(
  W,
  n_boot = 200L,
  strata = NULL,
  resample = c("stratified", "controls"),
  seed = NULL
)
```

## Arguments

- W:

  A matrix of treatment indicators with the values 0 and 1. Units are in
  rows.

- n_boot:

  The number of bootstrap draws.

- strata:

  An optional vector with one value for each unit. The values must not
  be `NA`.

- resample:

  Use `"stratified"` or `"controls"`. With `"controls"`, the treated
  units do not change in a draw. Refer to the section Inference of
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md).

- seed:

  An optional integer seed. The function does not change the state of
  the random number generator in the R session.

## Value

An integer matrix with `n_boot` rows and one column for each unit. Each
value is a row number of `W`.

## Details

A unit is a treated unit if its row of `W` has at least one value of 1.
The function draws the treated units and the other units separately.
With `strata`, it draws the units of each stratum separately. Thus, each
draw has the same number of treated units and of other units in each
stratum as the data.

A fit with fixed tuning parameters makes the same draws from the same
`seed`. Cross-validation with `cv_cells` uses random numbers before the
bootstrap. Then the draws are different.

A fit uses a treated unit only if the unit has a treated cell with an
outcome. If all treated cells of a unit have an `NA` outcome, the fit
draws the unit with the other units. Use `W` with such units set to 0 to
get the same draws.

## See also

[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)

## Examples

``` r
W <- matrix(0L, 8, 5)
W[1:3, 4:5] <- 1L
trop_boot_rows(W, n_boot = 3, seed = 1)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    1    3    1    5    8    6    5    6
#> [2,]    3    1    1    8    5    5    4    8
#> [3,]    1    1    1    8    8    5    5    4
trop_boot_rows(W, n_boot = 3, strata = rep(c("a", "b"), 4), seed = 1)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    1    3    2    5    7    4    8    8
#> [2,]    3    3    2    5    5    4    4    6
#> [3,]    3    3    2    5    5    8    4    4
```
