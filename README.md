# fasttrop

<!-- badges: start -->
[![R-CMD-check](https://github.com/farach/fasttrop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/farach/fasttrop/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

fasttrop is an R package for the Triply Robust Panel (TROP) estimator of
Athey, Imbens, Qu, and Viviano (2026). The estimator calculates the average
treatment effect on the treated (ATT) for a binary treatment in panel data.
All calculations run in C. The C code can use more than one thread with
OpenMP, and it uses the BLAS and LAPACK libraries of R. The authors of the
method did not write this package.

## How the estimator works

A cell is one unit in one period. For each treated cell, the estimator
calculates the outcome without the treatment. A model with two-way fixed
effects and a low-rank term calculates this outcome from the untreated cells. A
nuclear-norm penalty controls the low-rank term. Unit weights give more weight
to units with an outcome history that is similar to the treated unit. Time
weights give more weight to periods near the treated period.

The effect for a treated cell is the difference between the outcome of the
cell and the calculated outcome. The ATT is the mean of these effects.

## Reasons for this package

On 30 October 2026, CRAN had no package for the TROP estimator. These
implementations were available:

- The authors of the method wrote the Python package `trop`. It uses CVXPY.
- The Python package `diff-diff` has the local method and the global method.
- The R package `cfcompare` on GitHub has the local method in R code only.

The selection of the tuning parameters takes the most time. The leave-one-out
criterion calculates a new fit for each untreated cell at each candidate value.
The bootstrap calculates a new fit for each draw. On the test computer, the
default cross-validation of `cfcompare` took 17 minutes for a panel with 20
units and 15 periods. For a panel with 50 units and 40 periods, it took 27
minutes.

A research council selected the TROP estimator for this package. The council
used large language models as specialists, as a skeptic, and as judges. It
examined the R packages for difference-in-differences and the statistics and
econometrics papers from 2023 to 2026. The council accepted the package with
five conditions. The conditions are about agreement with other
implementations, tuning and inference, simulations, speed, and the limits of
the claims.

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

On the test computer with 12 threads, this example takes approximately 14
seconds. With fixed tuning parameters, the example takes approximately 5
seconds. For example, set `lambda = c(time = 0.1, unit = 0.5, nn = 1)`. Eq. 5
of the paper uses all 760 untreated cells. With all cells, the
cross-validation takes approximately 100 seconds.

## What the package contains

The package has two methods:

- The local method (`method = "local"`) is the estimator of Eq. 12 and
  Algorithm 2 of the paper. Each treated cell has its own unit weights, time
  weights, and fit. The leave-one-out criterion of Eq. 4 and Eq. 5 selects the
  tuning parameters. The search follows footnote 2 of the paper. You can use
  this method for block, staggered, and general designs.
- The global method (`method = "global"`) is the block estimator of the
  authors' package `trop`. It calculates one weighted fit. The time weights
  decrease with the distance from the center of the treated block. The unit
  weights decrease with the distance from the average treated unit before the
  treatment. Placebo cross-validation on untreated units selects the tuning
  parameters. You can use this method only for a block design.

In a block design, all treated units start the treatment in the same period
and stay treated. In a staggered design, the start periods are different, but
the treatment does not stop. In a general design, the treatment stops for one
or more units.

With `lambda = c(time = 0, unit = 0, nn = Inf)`, both methods give the same
result as imputation with two-way fixed effects (TWFE). The tests make sure of
this result.

The unit bootstrap of Algorithm 3 gives the standard error. The bootstrap
keeps the tuning parameters fixed. With `resample = "stratified"`, the
bootstrap draws treated and untreated units separately. With
`resample = "controls"`, the treated units do not change. The default
confidence interval is a normal interval. To get the percentile interval of
the authors' package, set `ci = "percentile"`.

Theorem 5.1 of the paper gives a bound for the bias of the estimator. The proof
is only for a block design. For other designs, fasttrop gives a warning
message. The paper also says that the bootstrap is correct only with many
treated units. Thus, fasttrop gives a warning message for a bootstrap with
fewer than five treated units.

## Validation

The research council set the conditions C1 to C5 before the build. This
section gives the results for C1, C2, and C3. The folder `inst/validation` has
the scripts, the test data, and all results.

### Agreement with the reference implementations (C1)

The tests compare fasttrop with two reference implementations at fixed
commits. The comparison uses four test panels and five sets of tuning
parameters. The table gives the largest absolute difference for each
comparison.

| Comparison | Largest difference |
|---|---|
| Local ATT, fasttrop and `diff-diff` | 4.0e-8 |
| Local effect for each cell, fasttrop and `diff-diff` | 1.1e-7 |
| Global ATT, fasttrop and `diff-diff` | 2.1e-8 |
| Global ATT, fasttrop and the authors' package | 2.5e-5 |
| Objective value with the same weights, relative to `diff-diff` | 7e-15 |

The difference from the authors' package is equal to the difference between
the two reference implementations. The cause is the lower accuracy of the
conic solver in that package. With the same weights, the objective values of
fasttrop and `diff-diff` are the same to machine precision. Thus, both solvers
find the same minimum.

### Tuning and inference (C2)

The tests make sure of these properties:

- The leave-one-out criterion does not use a cell in its own fit or in its own
  unit distances.
- Each bootstrap draw is equal to a new estimate on the resampled panel.
- The number of threads does not change the results.
- The local weights are equal to the weights from an independent R calculation
  of the formula in the paper.

The test suite has 472 expectations. On the test computer,
`R CMD check --as-cran --no-manual` gives the status OK.

### Simulations with a known effect (C3)

The script `inst/validation/known_truth.R` has three Monte Carlo simulations.
In each simulation, the true effect is 1 for each treated cell. The untreated
outcomes have unit effects, time effects, a factor with a linear trend, a
second factor, and noise. The treatment starts in the same period for all
treated units. The intervals come from the bootstrap of the package. In the
table, SE is the standard error and RMSE is the root mean square error.

| Simulation | Method | Bias (Monte Carlo SE) | RMSE | Coverage, normal interval | Coverage, percentile interval |
|---|---|---|---|---|---|
| A: 10 of 40 units, random selection | TWFE | -0.10 (0.05) | 0.75 | 0.890 | 0.870 |
| | TROP local | -0.04 (0.02) | 0.30 | 0.940 | 0.890 |
| | TROP global | -0.03 (0.02) | 0.29 | 0.925 | 0.930 |
| A2: 30 of 80 units, random selection | TWFE | -0.08 (0.03) | 0.31 | 0.920 | 0.930 |
| | TROP local | -0.04 (0.02) | 0.17 | 0.970 | 0.960 |
| | TROP global | -0.03 (0.02) | 0.16 | 0.990 | 0.980 |
| B: 10 of 40 units, selection on the trend | TWFE | 1.20 (0.05) | 1.37 | 0.500 | 0.500 |
| | TROP local | 0.64 (0.02) | 0.72 | 0.475 | 0.390 |
| | TROP global | 0.37 (0.02) | 0.47 | 0.700 | 0.645 |

With random selection, the bias of each estimator is less than three Monte
Carlo standard errors. The RMSE of TROP is 44% to 62% less than the RMSE of
TWFE. All three methods use the same panels. Thus, the small negative biases
are the same noise for all methods. A separate run with 2,000 replications of
TWFE gave a bias of 0.022 for simulation A and -0.001 for simulation A2. The
Monte Carlo standard errors were 0.016 and 0.007.

For TROP, the coverage of the normal interval is near 95%. With 100 draws, the
coverage of the percentile interval is only 89% for the local method in
simulation A. For this reason, the normal interval is the default. In
simulation A2 with 30 treated units, no coverage is less than 95% by more than
two Monte Carlo standard errors. The global interval has a coverage of 0.99,
so it is wider than necessary. The paper also says that the bootstrap is
correct only with many treated units.

In simulation B, the selection of treated units uses their loading on the
trend factor. Thus, TWFE has a large bias. TROP removes approximately half of
this bias with the local method and two thirds with the global method. But
some bias stays, and the coverage of the intervals is less than 95%.

This bias does not come from the solver. With the default tolerance and with a
tolerance of 1e-10, the estimates have a difference of less than 6e-4. Pure
matrix completion with a very small penalty has the same bias. In this
simulation, the treated units have the largest trends. Thus, no weights on
untreated units or on earlier periods can make a good match. TROP cannot
remove a bias from trends that the untreated units do not have.

## Speed

The request for this package was speed equal to C code. All calculations are
in the file `src/trop_core.c`, which does not use the R API. The R code
examines the inputs, makes the random draws, and makes the result object.

To measure the cost of the R code, the C program `inst/bench/trop_driver.c`
uses the same C core. The program uses the same compiler, flags, BLAS, and
LAPACK as R. It runs the same plan as R, with the same sampled cells, placebo
sets, and bootstrap rows. Its results are the same as the results of the
package, bit for bit.

| Workload | Threads | C program (s) | `trop_run()` (s) | `trop_matrix()` (s) |
|---|---|---|---|---|
| Local, fixed lambda, no bootstrap, 20 x 12 | 1 | 0.0204 | 0.0207 | 0.0221 |
| | 4 | 0.0056 | 0.0060 | 0.0064 |
| Local, fixed lambda, 100 bootstrap draws, 30 x 20 | 1 | 14.38 | 14.10 | 15.70 |
| | 4 | 3.99 | 3.88 | 3.75 |
| Local, CV on 60 cells, 50 draws, 30 x 20 | 1 | 4.80 | 4.34 | 4.69 |
| | 4 | 1.12 | 1.17 | 1.15 |
| Global, placebo CV with 50 sets, 200 draws, 50 x 30 | 1 | 7.03 | 7.21 | 7.81 |
| | 4 | 1.86 | 1.85 | 1.89 |

The values are medians of three rounds. In each round, the three paths ran in
a random sequence. The function `trop_run()` starts the C calculations for a
prepared plan. The function `trop_matrix()` is the full function call, with
the input checks and the random draws.

For the paired runs of one second or more, the ratio of R time to C time had
these values:

- For `trop_run()`, the geometric mean was 0.996, with a 95% interval of 0.970
  to 1.023.
- For `trop_matrix()`, the geometric mean was 1.012, with a 95% interval of
  0.983 to 1.041.

On this laptop, single runs had a variation of approximately 10%. In the
limits of this variation, the R package has the same speed as the C program.
The fixed cost of the R code is approximately 0.3 ms for each call. This cost
is important only for very small panels. The C program and R gave the same
tuning parameters, estimates, and bootstrap draws. Runs with 1 thread and with
4 threads also gave the same results.

### Windows power throttling

The first version of this benchmark showed that R was approximately 40% slower
than the C program. The cause was not in fasttrop. The same LAPACK calls were
1.34 times slower in a long R session than in a new process. R was slower in 28
of 30 pairs with a random sequence.

Windows uses power throttling (EcoQoS) on processes that run in the background
for a long time. The slowdown increased with the time that the R process was
busy. The slowdown stopped when the process was not under power throttling.
The short C program stopped before Windows started power throttling.

fasttrop now tells Windows not to use power throttling on its threads while
the C code runs. When the C code stops, Windows controls the threads again. To
keep the usual Windows policy, use `trop_control(throttle = TRUE)`. On other
operating systems, this option has no effect.

Two limits apply to this comparison:

- The comparison uses the same C core. A different algorithm or a faster BLAS,
  for example OpenBLAS, can be faster. When R uses a faster BLAS, fasttrop also
  uses it.
- The benchmarks ran on a Snapdragon X Elite laptop with x64 R 4.6.1 under
  emulation on Windows on ARM. The package and the C program both ran under
  emulation. Thus, the conditions of the comparison were the same. A native
  build can be faster.

### Comparison with the other implementations

The other implementations used the same two panels and the same fixed tuning
parameters, (time, unit, nn) = (0.3, 0.8, 0.5). Each package used its default
solver options. The scripts are `inst/bench/reference_timings.R` and
`inst/bench/reference_timings.py`.

| Task | Panel | fasttrop, 1 thread | fasttrop, 4 threads | Authors' `trop` | `diff-diff` | `cfcompare` |
|---|---|---|---|---|---|---|
| Local ATT | 20 x 15 | 0.044 s | 0.016 s | | 3.42 s | 5.93 s |
| Local ATT | 40 x 20 | 0.33 s | 0.088 s | | 9.33 s | 49.7 s |
| Global ATT | 20 x 15 | 0.0032 s | 0.0040 s | 0.100 s | 0.551 s | |
| Global ATT | 40 x 20 | 0.010 s | 0.0084 s | 0.252 s | 0.951 s | |
| One evaluation of the leave-one-out criterion, 284 cells | 20 x 15 | 0.85 s | 0.20 s | | 56.6 s | |

With one thread, fasttrop was 25 to 31 times faster than the authors' package
for the global method. It was 29 to 174 times faster than `diff-diff` and 136
to 152 times faster than `cfcompare`. The ATTs of fasttrop and the two Python
packages had a difference of less than 1e-5 at their default tolerances. The
ATTs of `cfcompare` had a difference of up to 0.006. Thus, `cfcompare` does
not calculate exactly the same value.

`diff-diff` ran in Python only, because its optional Rust code was not
available. With the Rust code, its times can be shorter.

The cross-validation procedures of the packages are different, but these
times give context. For a different panel with 20 units and 15 periods, the
default cross-validation of `cfcompare` took 1,023 seconds. It used 120
sampled cells and its own grid. On the same panel, the default
cross-validation of fasttrop used all 284 untreated cells. It took 31 seconds
with one thread and 3.8 seconds with 12 threads.

## Limits

These limits apply to version 0.1.0:

- The package cannot use covariates.
- Outcomes can be `NA`, and the fits do not use these cells. The tests for the
  tuning and for the bootstrap use only a small number of `NA` cells.
- The equations come from arXiv version 3 of the paper. A line-by-line
  comparison with the journal version is not complete.
- GitHub Actions runs `R CMD check` on Linux, macOS, and Windows. The
  validation scripts and the benchmark scripts ran only on Windows.

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
