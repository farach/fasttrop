# Tuning grids

`trop_grid()` gives the candidate values of the three tuning parameters.
These parameters are the decay of the time weights, the decay of the
unit weights, and the nuclear-norm penalty. A value of 0 for `time` or
`unit` gives equal weights. A value of `Inf` for `nn` removes the
low-rank term. Then the model has only two-way fixed effects.

## Usage

``` r
trop_grid(
  time = c(0, 0.1, 0.5, 1, 2, 5),
  unit = c(0, 0.1, 0.5, 1, 2, 5),
  nn = c(0.01, 0.1, 1, 10, Inf)
)
```

## Arguments

- time:

  The candidate values for the decay of the time weights. The values
  must be finite and 0 or more.

- unit:

  The candidate values for the decay of the unit weights. The values
  must be finite and 0 or more.

- nn:

  The candidate values for the nuclear-norm penalty. The values must be
  0 or more. You can use `Inf`.

## Value

A list of class `trop_grid` with the numeric vectors `time`, `unit`, and
`nn`.

## Examples

``` r
trop_grid(nn = c(0.1, 1, Inf))
#> $time
#> [1] 0.0 0.1 0.5 1.0 2.0 5.0
#> 
#> $unit
#> [1] 0.0 0.1 0.5 1.0 2.0 5.0
#> 
#> $nn
#> [1] 0.1 1.0 Inf
#> 
#> attr(,"class")
#> [1] "trop_grid"
```
