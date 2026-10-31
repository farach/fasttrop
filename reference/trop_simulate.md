# Simulated panel with a factor structure

`trop_simulate()` makes a balanced panel with unit effects, time
effects, a low-rank interactive term, and Gaussian noise. Each unit gets
a score, which is `selection` times its first factor loading plus random
noise. The units with the highest scores get the treatment.

With a positive value of `selection`, a unit with a larger loading has a
higher probability of treatment. With a negative value, a unit with a
smaller loading has a higher probability. With zero, the loading has no
effect on the selection. Imputation with two-way fixed effects (TWFE
imputation) does not model this factor, so its estimate for one panel
can have a large error. TROP can reduce such errors, but some error can
stay.

## Usage

``` r
trop_simulate(
  n_units = 40L,
  n_periods = 30L,
  n_treated = 8L,
  n_post = 5L,
  effect = 1,
  rank = 2L,
  sigma = 0.5,
  selection = 1,
  design = c("block", "staggered"),
  seed = NULL
)
```

## Arguments

- n_units, n_periods:

  The number of units and the number of periods.

- n_treated:

  The number of treated units.

- n_post:

  For a block design, the number of treated periods at the end of the
  panel. For a staggered design, the treatment starts in one of the last
  `n_post` periods.

- effect:

  The constant treatment effect on each treated cell.

- rank:

  The rank of the interactive term.

- sigma:

  The standard deviation of the noise.

- selection:

  The weight of the first factor loading in the selection of the treated
  units. With zero, the selection is random.

- design:

  Use `"block"` or `"staggered"`.

- seed:

  An optional seed. The function does not change the state of the random
  number generator in the R session.

## Value

A data frame with the columns `unit`, `time`, `y`, `d`, and `y0`. The
column `y0` contains the outcome without the treatment.

## Examples

``` r
head(trop_simulate(seed = 1))
#>   unit time          y d         y0
#> 1    1    1  0.2388666 0  0.2388666
#> 2    2    1  0.5934588 0  0.5934588
#> 3    3    1 -2.5737706 0 -2.5737706
#> 4    4    1  1.2415187 0  1.2415187
#> 5    5    1  0.4650530 0  0.4650530
#> 6    6    1  0.5049260 0  0.5049260
```
