# fasttrop

fasttrop estimates the average effect of a binary treatment in panel data with
the Triply Robust Panel (TROP) estimator of Athey, Imbens, Qu, and Viviano
(2026). The estimator fits two-way fixed effects plus a nuclear-norm penalized
low-rank term on untreated cells, with unit weights based on how closely each
unit's untreated history tracks the treated unit and time weights that decay
with distance in periods. The numerical work runs in C with OpenMP and calls
the BLAS and LAPACK that R itself uses.

This is an independent implementation. It is not affiliated with the authors
of the method.

## Why this package

On 30 October 2026 no CRAN package implemented TROP. The authors' reference
code is a Python package built on CVXPY, `diff-diff` implements it in Python,
and `cfcompare` implements it in pure R on GitHub. Tuning is the expensive
step. The leave-one-out criterion refits a weighted low-rank model for every
untreated cell at every candidate tuning value, and the bootstrap repeats the
fit for every draw. On the test machine, `cfcompare` needed 17 minutes for its
default cross-validation on a 20 x 15 panel and 27 minutes on a 50 x 40 panel.

The method was chosen by a multi-model research council that reviewed the R
packages for difference-in-differences and the 2023 to 2026 statistics and
econometrics literature. The council approved the build under five conditions
covering agreement with reference implementations, tuning and inference,
simulations with a known answer, speed, and the scope of claims. The
validation and speed sections below report the results.

## Installation

fasttrop is not on CRAN. Install it from the package folder with a C compiler
available (Rtools on Windows):

```r
install.packages("path/to/fasttrop", repos = NULL, type = "source")
```

## Example

```r
library(fasttrop)
panel <- trop_simulate(n_units = 40, n_periods = 20, n_treated = 8, n_post = 5, seed = 1)

# Tuning by leave-one-out cross-validation on 200 sampled untreated cells,
# then 200 bootstrap draws.
fit <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
            method = "local", cv_cells = 200, n_boot = 200, seed = 1)
summary(fit)
confint(fit)
head(as.data.frame(fit)) # effects for each treated cell
```

On the test machine with 12 threads this example takes about 14 seconds. With
fixed tuning parameters, for example `lambda = c(time = 0.1, unit = 0.5, nn = 1)`,
it takes about 5 seconds. Cross-validation over all 760 untreated cells, as in
the paper's Eq. 5, takes about 100 seconds.

## What is implemented

`method = "local"` is the estimator in the paper (Eq. 12 and Algorithm 2). Each
treated cell gets its own unit and time weights and its own fit. The tuning
parameters minimize the leave-one-out criterion of Eq. 4 and 5, found with the
coordinate search described in the paper's footnote 2. It works for block,
staggered, and general assignment. Setting `cv_cells` evaluates the criterion
on a random sample of untreated cells instead of all of them.

`method = "global"` is the block estimator in the authors' Python package. It
makes one weighted fit, with time weights centred on the treated block and unit
weights based on distance to the average treated unit before treatment. Tuning
uses placebo cross-validation on untreated units.

Standard errors come from the unit bootstrap of Algorithm 3 with the tuning
parameters held fixed. Intervals are normal by default, and `ci = "percentile"`
gives the percentile interval used by the authors' package.
`resample = "stratified"` draws treated and untreated units separately, and
`resample = "controls"` keeps the treated units fixed. With
`lambda = c(time = 0, unit = 0, nn = Inf)` the estimator reduces to two-way
fixed effects imputation.

The paper proves its bias bound (Theorem 5.1) for block assignment only.
fasttrop warns when the design is staggered or when treatment switches off, and
when fewer than five units are treated, because the bootstrap's validity
argument needs a growing number of treated units.

## Validation

Agreement with the references (condition C1) was checked on four fixture panels
and five tuning settings against two independent implementations pinned to
fixed commits. The scripts, fixtures, and detailed results are in
`inst/validation`. Local effects match `diff-diff` to within 4e-8 for the ATT and
1.1e-7 for single cells. Global effects match `diff-diff` to within 2.1e-8 and
the authors' package to within 2.5e-5, which is the accuracy of the authors'
conic solver. Given identical weights, fasttrop's objective value equals
`diff-diff`'s to machine precision.

The tests for condition C2 confirm that the leave-one-out criterion leaves each
held-out cell out of its own fit and its own unit distances, that every
bootstrap draw equals a fresh estimate on the resampled panel, and that results
do not depend on the number of threads. The local weights are also checked
against an independent R implementation of the paper's definition. The test
suite has 467 expectations and passes `R CMD check --as-cran` with no notes.

Condition C3 asked for simulations with a known answer, under block
assignment and a factor model. `inst/validation/known_truth.R` runs three
designs with a constant effect of 1. Untreated outcomes have unit and time
effects, a linear-trend factor, a second factor, and noise. Each design
compares TROP with two-way fixed effects (TWFE), and all intervals come from
the package's own bootstrap.

| Design | Method | Bias (MC SE) | RMSE | Normal CI coverage | Percentile CI coverage |
|---|---|---|---|---|---|
| A: 10 of 40 units treated at random | TWFE | -0.10 (0.05) | 0.75 | 0.890 | 0.870 |
| | TROP local | -0.04 (0.02) | 0.30 | 0.940 | 0.890 |
| | TROP global | -0.03 (0.02) | 0.29 | 0.925 | 0.930 |
| A2: 30 of 80 units treated at random | TWFE | -0.08 (0.03) | 0.31 | 0.920 | 0.930 |
| | TROP local | -0.04 (0.02) | 0.17 | 0.970 | 0.960 |
| | TROP global | -0.03 (0.02) | 0.16 | 0.990 | 0.980 |
| B: 10 of 40 treated, selected on the trend loading | TWFE | 1.20 (0.05) | 1.37 | 0.500 | 0.500 |
| | TROP local | 0.64 (0.02) | 0.72 | 0.475 | 0.390 |
| | TROP global | 0.37 (0.02) | 0.47 | 0.700 | 0.645 |

With random assignment, the three estimators agree with the true effect to
within two to three Monte Carlo standard errors, and TROP's RMSE is 44% to 62%
lower than TWFE's. The small negative biases are sampling noise shared by all
three methods, which see the same panels. A separate run of 2,000 TWFE
replications gave biases of 0.022 (0.016) for design A and -0.001 (0.007) for
design A2. Normal intervals built from the bootstrap standard error cover close
to 95% for TROP. Percentile intervals from 100 draws cover only 89% for the
local method in design A, which is why `ci = "normal"` is the default. With 30
treated units instead of 10 (design A2), coverage is at or above 95% within
Monte Carlo error for every method, and the global method's intervals are
conservative (0.99). This fits the paper's requirement of a growing number of
treated units.

When treatment is selected on the loading of the trending factor (design B),
TWFE is badly biased. TROP removes about half of that bias with the local
method and two thirds with the global method, but not all of it, and the
intervals undercover because they ignore the remaining bias. This is not a
solver artifact. Estimates at the default tolerance and at 1e-10 differ by less
than 6e-4, and pure matrix completion with a tiny penalty shows the same bias.
In this design the treated units have the steepest trends, so no weighting of
untreated units or earlier periods can balance them, and the bound in
Theorem 5.1 is not small. Users should not expect TROP to remove confounding
that the untreated units cannot reproduce.

## Speed, and whether it matches C

The request behind this package was to make the code run as fast as its C
equivalent. All numerical work, including the cross-validation loops and the
bootstrap, runs in `src/trop_core.c`, which uses no R API. R validates inputs,
draws the random numbers, and assembles the result. To test whether the R layer
costs anything, `inst/bench/trop_driver.c` compiles the same core into a
standalone program with R's compiler, flags, BLAS, and LAPACK. It runs the
same plan as R (the same sampled cells, placebo sets, and bootstrap rows), and
its results are bit-identical to the package's.

| Workload | Threads | C program (s) | R `trop_run()` (s) | R `trop_matrix()` (s) |
|---|---|---|---|---|
| Local, fixed lambda, no bootstrap, 20 x 12 | 1 | 0.0204 | 0.0207 | 0.0221 |
| | 4 | 0.0056 | 0.0060 | 0.0064 |
| Local, fixed lambda, 100 bootstrap draws, 30 x 20 | 1 | 14.38 | 14.10 | 15.70 |
| | 4 | 3.99 | 3.88 | 3.75 |
| Local, CV on 60 cells, 50 draws, 30 x 20 | 1 | 4.80 | 4.34 | 4.69 |
| | 4 | 1.12 | 1.17 | 1.15 |
| Global, placebo CV with 50 sets, 200 draws, 50 x 30 | 1 | 7.03 | 7.21 | 7.81 |
| | 4 | 1.86 | 1.85 | 1.89 |

Values are medians of three rounds, with the three paths run in random order
within each round. `trop_run()` makes the native calls for a prepared plan, and
`trop_matrix()` is the full call including input checks and random draws.
Across the 18 paired runs that took at least a second, the time of
`trop_run()` divided by that of the C program had a geometric mean of 0.996
(95% interval 0.970 to 1.023), and the same ratio for `trop_matrix()` was 1.012
(0.983 to 1.041). Single runs on this laptop varied by about 10%, and within
that noise the R package runs at the speed of the C program. The fixed cost of
the R layer shows only on the smallest problem, at about 0.3 ms per call. The C
program and R return identical tuning parameters, estimates, and bootstrap
draws, and so do runs with 1 and 4 threads.

The first version of this benchmark found the R path about 40% slower than the
C program. The cause was not in fasttrop. Identical LAPACK calls ran a median
1.34 times slower inside a long-running R session than in a fresh process
(R slower in 28 of 30 randomized pairs). The slowdown grew the longer the R
process stayed busy, and it disappeared when the process opted out of Windows
power throttling (EcoQoS), which Windows applies to long-running background
processes. The short-lived C program finished before it was throttled.
fasttrop now asks Windows not to throttle the threads that run its kernels
while they run, and hands the decision back to the system afterwards.
`trop_control(throttle = TRUE)` keeps the system policy, and the option does
nothing on other platforms.

Two limits apply to the parity claim. Parity is measured against the same C
core, so a different algorithm or an optimized BLAS such as OpenBLAS would be
faster, and fasttrop gains the same speedup when R is linked to that BLAS. The
benchmarks also ran on a Snapdragon X Elite laptop with x64 R 4.6.1 under
Windows on ARM emulation. Both the package and the C program were emulated, so
the comparison is fair, but native builds would be faster in absolute terms.

The other implementations were timed on the same two panels with the same
fixed tuning parameters, (time, unit, nn) = (0.3, 0.8, 0.5), and each package's
default solver settings (`inst/bench/reference_timings.R` and `.py`).

| Task | Panel | fasttrop, 1 thread | fasttrop, 4 threads | Authors' `trop` | `diff-diff` | `cfcompare` |
|---|---|---|---|---|---|---|
| Local ATT | 20 x 15 | 0.044 s | 0.016 s | | 3.42 s | 5.93 s |
| Local ATT | 40 x 20 | 0.33 s | 0.088 s | | 9.33 s | 49.7 s |
| Global ATT | 20 x 15 | 0.0032 s | 0.0040 s | 0.100 s | 0.551 s | |
| Global ATT | 40 x 20 | 0.010 s | 0.0084 s | 0.252 s | 0.951 s | |
| One leave-one-out evaluation, 284 cells | 20 x 15 | 0.85 s | 0.20 s | | 56.6 s | |

On one thread fasttrop was 25 to 31 times faster than the authors' package for
the global estimator, 29 to 174 times faster than `diff-diff`, and 136 to 152
times faster than `cfcompare`. Its estimates agreed with the two Python packages
to within 1e-5 at their default tolerances. `cfcompare`'s estimates differed by
up to 0.006, so it computes nearly but not exactly the same quantity.
`diff-diff` ran as pure Python because its optional Rust backend was not
available here, so these numbers understate that package's speed. For context
only, since the procedures differ: `cfcompare`'s default cross-validation (120
sampled cells and its own grid) took 1,023 seconds on a 20 x 15 panel, while
fasttrop's default cross-validation over all 284 untreated cells of the same
panel took 31 seconds on one thread and 3.8 seconds on 12 threads.

## Limitations

- Covariates are not supported.
- Missing outcomes are allowed and excluded from the fits, but the tuning and
  bootstrap code has only been tested with a few missing cells.
- The equations were checked against arXiv version 3 of the paper. The journal
  version of record has not been compared line by line.
- The package has been built and tested on Windows only.

## References

Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026). Triply Robust Panel
Estimators. *Journal of Applied Econometrics*. doi:10.1002/jae.70061.
arXiv:2508.21536.

Reference implementations used for validation: `trop` by the paper's authors
(github.com/zhaonanq/TROP-Estimator, MIT), `diff-diff`
(github.com/igerber/diff-diff, MIT), and `cfcompare` by Takuma Iwasaki
(github.com/takuma1102/cfcompare, MIT).
