# fasttrop <img src="man/figures/logo.png" align="right" height="139" alt="fasttrop hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/farach/fasttrop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/farach/fasttrop/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

fasttrop is an R package for the Triply Robust Panel (TROP) estimator of
Athey, Imbens, Qu, and Viviano (2026). The estimator gives the effect of a
binary treatment in panel data. An example is the effect of a new law on the
states that adopt it. The package calculates the average treatment effect on
the treated (ATT), a standard error, and a confidence interval. The model fits
run in C, and they can use more than one thread.

The authors of the method did not write this package.

The [fasttrop website](https://farach.github.io/fasttrop/) explains the method
in plain words. It also shows how to read the results, the limits of the
method, and how we checked the package.

## What TROP does

Panel data has many units in many periods. The units can be states, firms, or
people. A cell is one unit in one period. A treated cell is a cell with the
treatment.

For each treated cell, TROP calculates the outcome of the cell without the
treatment. It uses only the untreated cells for this calculation. The effect
for the cell is the difference between the observed outcome and the calculated
outcome. The ATT is the mean of these effects.

TROP combines three ideas:

- Unit weights give more weight to units with an outcome history that is
  similar to the treated unit.
- Time weights give more weight to periods near the treated period.
- A model finds a small number of patterns that the units share over time. A
  penalty controls how many patterns the model uses.

Cross-validation selects the strength of each idea from the data. Under the
assumptions of the paper, the bias of TROP has a limit. This limit is the
product of three terms: the unit imbalance, the time imbalance, and the bias
of the model. If one term goes to zero and the other terms stay limited, the
limit goes to zero. This is the meaning of "triply robust". The proof is only
for a block design with fixed weights.

## Why fasttrop

We wanted a TROP estimator for R that is fast and that agrees with the
reference implementations. The tuning step and the bootstrap calculate
thousands of fits. Thus, speed is important.

When we started fasttrop, these implementations were available:

- The authors of the method wrote the Python package `trop`. It has the global
  method for a block design.
- The Python package `diff-diff` has the local method and the global method.
- The R package `cfcompare` on GitHub has the local method in R code.

fasttrop has the local method and the global method. In our benchmarks with
fixed tuning parameters and one thread, fasttrop was 25 to 174 times faster
than the Python packages. It was 136 to 152 times faster than `cfcompare`.

## Install the package

fasttrop is not on CRAN at this time. To install the package from GitHub, do
these steps:

1. Make sure that you have a C compiler. On Windows, Rtools has the compiler.
2. If you do not have the `remotes` package, install it.
3. Run this R code:

```r
remotes::install_github("farach/fasttrop")
```

## Example

```r
library(fasttrop)
panel <- trop_simulate(n_units = 40, n_periods = 20, n_treated = 8, n_post = 5, seed = 1)

# Select the tuning parameters with leave-one-out cross-validation on 200
# sampled untreated cells. Then calculate 200 bootstrap draws.
fit <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
            method = "local", cv_cells = 200, n_boot = 200, seed = 1)
summary(fit)
confint(fit)
head(as.data.frame(fit)) # The effect for each treated cell
```

The cross-validation takes most of the time. To make the example faster, set
fixed tuning parameters, for example `lambda = c(time = 0.1, unit = 0.5, nn = 1)`.
With `cv_cells = NULL`, the cross-validation uses all untreated cells, as in
Eq. 5 of the paper. This takes more time.

## Methods

- The local method (`method = "local"`) is the estimator of Eq. 12 and
  Algorithm 2 of the paper. Each treated cell has its own unit weights, time
  weights, and fit. The leave-one-out criterion of Eq. 4 and Eq. 5 selects the
  tuning parameters. You can use this method for block, staggered, and general
  designs.
- The global method (`method = "global"`) is the block estimator of the
  authors' package `trop`. It calculates one weighted fit for all treated
  cells. Placebo cross-validation on untreated units selects the tuning
  parameters. You can use this method only for a block design.

In a block design, all treated units start the treatment in the same period
and stay treated. In a staggered design, the start periods are different, but
the treatment does not stop. In a general design, the treatment stops for one
or more units.

A unit bootstrap gives the standard error. The resampling follows Algorithm 3
of the paper. fasttrop keeps the selected tuning parameters fixed in each
draw. The default confidence interval is a normal interval. To get the
percentile interval of the authors' package, set `ci = "percentile"`.

With `lambda = c(time = 0, unit = 0, nn = Inf)`, both methods give the same
result as imputation with two-way fixed effects (TWFE). The tests make sure of
this result.

## Validation

We compared fasttrop with two reference implementations at fixed commits. The
comparison used four test panels and five sets of tuning parameters. The table
gives the largest absolute difference for each comparison.

| Comparison | Largest difference |
|---|---|
| Local ATT, fasttrop and `diff-diff` | 4.0e-8 |
| Local effect for each cell, fasttrop and `diff-diff` | 1.1e-7 |
| Global ATT, fasttrop and `diff-diff` | 2.1e-8 |
| Global ATT, fasttrop and the authors' package | 2.5e-5 |

The difference from the authors' package is equal to the difference between
the two reference implementations. The cause is the lower accuracy of the
conic solver in that package.

The tests also examine the tuning and the bootstrap. For example, the
leave-one-out criterion does not use a cell in its own fit, and the number of
threads does not change the results. The test suite has 472 expectations.

Monte Carlo simulations with a known effect gave these results:

- With random selection of the treated units, the root mean square error
  (RMSE) of TROP was 44% to 62% less than the RMSE of TWFE. The coverage of
  the normal intervals of TROP was from 92.5% to 99%.
- When the selection of the treated units used their trend, TWFE had a large
  bias. TROP removed approximately half of this bias with the local method and
  two thirds with the global method. But some bias stayed, and the coverage of
  the intervals was less than 95%.

The folder `inst/validation` has the scripts, the data, and the saved results
of the main checks.

## Speed

We wanted fasttrop to be fast. The model fits are in the file
`src/trop_core.c`, which does not use the R API. The R code examines the
inputs and makes the random draws. It also calculates the standard error and
the interval from the bootstrap draws, and it makes the result object.

A benchmark compared the package with a C program that uses the same C core.
In the limits of the measurement, the R package had the same speed as the C
program. For a 20 by 12 panel with fixed tuning parameters and no bootstrap,
the full R call added 0.8 to 1.7 ms in our benchmark. A call with tuning or a
bootstrap can add more R time.

The table compares fasttrop with the other implementations. All runs used the
same panels and the same fixed tuning parameters. An empty cell shows that the
package does not have the method. Times depend on the computer.

| Task | Panel | fasttrop, 1 thread | Authors' `trop` | `diff-diff` | `cfcompare` |
|---|---|---|---|---|---|
| Local ATT | 20 x 15 | 0.044 s | | 3.42 s | 5.93 s |
| Local ATT | 40 x 20 | 0.33 s | | 9.33 s | 49.7 s |
| Global ATT | 20 x 15 | 0.0032 s | 0.100 s | 0.551 s | |
| Global ATT | 40 x 20 | 0.010 s | 0.252 s | 0.951 s | |

`diff-diff` ran in Python only, because its optional Rust code was not
available. With the Rust code, its times can be shorter. The ATTs of
`cfcompare` differed by up to 0.006, so it does not calculate exactly the same
value. On Windows, fasttrop tells Windows not to make its threads slower while
the C code runs. Refer to `?trop_control`. The folder `inst/bench` has the
scripts and the saved results of the benchmarks.

## Limits

These limits apply to version 0.1.0:

- Theorem 5.1 of the paper (triple robustness) has a proof only for a block
  design. For other designs, fasttrop gives a warning message.
- The bootstrap is correct only with many treated units. fasttrop gives a
  warning message for a bootstrap with fewer than five treated units.
- The global method works only for a block design.
- The package cannot use covariates.
- In our simulation in which the treated units had larger trends, TROP
  reduced the bias, but some bias stayed. This is the result of one
  simulation, not a general limit of the method.
- Outcomes can be `NA`, and the fits do not use these cells. The tests for the
  tuning and for the bootstrap use only a small number of `NA` cells.
- The equations come from arXiv version 3 of the paper. A line-by-line
  comparison with the journal version is not complete.
- We ran the validation scripts and the benchmark scripts only on Windows.
  GitHub Actions runs `R CMD check` on Linux, macOS, and Windows.

## Cite the method and the package

If you use fasttrop, cite the paper. To see the citations for the paper and
the package, run `citation("fasttrop")`.

## References

Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026). Triply Robust Panel
Estimators. *Journal of Applied Econometrics*. doi:10.1002/jae.70061.
arXiv:2508.21536.

The validation used these reference implementations:

- Python package `trop` from the authors of the paper,
  github.com/zhaonanq/TROP-Estimator (MIT license)
- Python package `diff-diff`, github.com/igerber/diff-diff (MIT license)
- R package `cfcompare` by Takuma Iwasaki, github.com/takuma1102/cfcompare
  (MIT license)

## License

The package has the MIT license. Refer to the file `LICENSE.md`.
