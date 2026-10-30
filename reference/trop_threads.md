# Default number of threads

`trop_threads()` gives the default number of threads for
[`trop_control()`](https://farach.github.io/fasttrop/reference/trop_control.md).

## Usage

``` r
trop_threads()
```

## Value

An integer. If the option `fasttrop.threads` has a value, the function
returns this value. If not, the function returns the number of physical
cores, with a minimum of 1. When the environment variable
`_R_CHECK_LIMIT_CORES_` has a value other than `"false"`, the maximum is
2. `R CMD check --as-cran` sets this variable.

## Examples

``` r
trop_threads()
#> [1] 2
```
