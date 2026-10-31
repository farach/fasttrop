# Validation

This folder has the scripts and the results for the validation of fasttrop.
The validation has five parts: agreement with the reference implementations,
tuning and inference, simulations with a known effect, speed, and the limits
of the claims. The package README gives a summary.

## Agreement with the reference implementations

The validation uses two reference implementations at fixed commits:

- The Python package `trop` from the authors of the paper
  (zhaonanq/TROP-Estimator, MIT license), commit
  `43738cdaa8eb8356566903bb010bbf6ce4cb31b9`. It has the global block
  estimator. It uses CVXPY with the CLARABEL solver.
- The Python package `diff-diff` (igerber/diff-diff, MIT license), commit
  `d46a5fd8b377923d944dc47bcf6d82233d0ebcfd`. It has the local estimator of
  Eq. 12 and Algorithm 2 of the paper, and the global estimator. It runs in
  Python only, with a tolerance of 1e-9.

The folder `tests/testthat/fixtures` has four simulated panels:

- D1: 12 units and 10 periods, block design
- D2: 20 units and 15 periods, block design
- D3: 15 units and 12 periods, staggered design
- D4: 16 units and 12 periods, block design with one treated unit

Each panel has five sets of tuning parameters. The sets are TWFE imputation,
weighted TWFE imputation, matrix completion, and two sets for the full TROP
estimator. The script
`data-raw/oracle_fixtures.py` calculates the reference results. fasttrop uses a
tolerance of 1e-10.

The script `oracle_report.R` calculates the largest absolute differences. The
files `oracle-att-comparison.csv` and `oracle-objective-comparison.csv` have
all values.

| Comparison | Largest difference |
|---|---|
| Local ATT and `diff-diff` | 4.0e-8 |
| Local effects for each cell and `diff-diff` | 1.1e-7 |
| Global ATT and `diff-diff` | 2.1e-8 |
| Global effects for each cell and `diff-diff` | 6.0e-8 |
| Global ATT and the authors' package `trop` | 2.5e-5 |
| Objective value with the same weights, relative to `diff-diff` | 7e-15 |
| Fitted values with the same weights and `diff-diff` | 2.5e-8 |

The difference from the authors' package is equal to the difference between
the two reference implementations. The cause is the lower accuracy of the
conic solver in that package. Without the low-rank term, the difference is
3e-12 or less.

With the same weights, the objective values of fasttrop and `diff-diff` are
the same to machine precision. Thus, both solvers find the same minimum. The
file `tests/testthat/test-oracle.R` makes sure of these results. The
tolerances are 1e-6 for `diff-diff`, 1e-4 for the authors' package, and a
relative 1e-7 for the objective value.

## Tuning and inference

The file `tests/testthat/test-tuning-inference.R` makes sure of these
properties:

- The leave-one-out criterion is equal to a direct calculation in which the
  treatment indicator of each cell of the criterion is 1. Thus, the criterion
  does not use a cell in its own fit or in its own unit distances.
- The placebo criterion is the RMSE of the placebo estimates.
- The coordinate searches give the score of the selected point. The global
  search stops at a fixed point, where no single move on the grid gives a
  better score.
- Each bootstrap draw is equal to a new estimate on the resampled panel, for
  both methods. fasttrop keeps the tuning parameters fixed in each draw.
  Algorithm 3 of the paper gives the resampling, but it does not say whether
  to select the tuning parameters again.
- The stratified bootstrap draws treated and untreated units separately. The
  option `resample = "controls"` keeps the treated units fixed.
- With a seed, the results are the same in each run. The functions do not
  change the state of the random number generator of the caller. The number
  of threads does not change the results.

The file `tests/testthat/test-estimator.R` compares the local weights with an
independent R calculation of the formula in the paper. Another test makes sure
that `lambda = c(time = 0, unit = 0, nn = Inf)` gives the same result as
TWFE imputation. This test has cells with `NA` outcomes.

Two more tests in this file compare TWFE imputation with a conventional TWFE
regression. In a balanced block design, the two results are the same. In a
staggered design with no noise, TWFE imputation gives the true ATT, and the
regression does not.

## Simulations with a known effect in a block design

The script `known_truth.R` runs the simulations A, A2, and B. The top of the
script gives the details of each simulation. Simulations A and B used 200
replications, and simulation A2 used 100 replications. The script
`summarize_known_truth.R` makes the files `known-truth-summary-*.csv` from the
saved draws. In these files, the method `twfe` is TWFE imputation. In a block
design with no `NA` outcomes, it gives the same estimate as a conventional
TWFE regression.

In each simulation, the true effect is 1 for each treated cell. The untreated
outcomes have unit effects, time effects, a factor with a linear trend, a
second factor, and noise. The treatment starts in the same period for all
treated units. The intervals come from the bootstrap of the package. In the
table, SE is the standard error and RMSE is the root mean square error.

| Simulation | Method | Bias (Monte Carlo SE) | RMSE | Coverage, normal interval | Coverage, percentile interval |
|---|---|---|---|---|---|
| A: 10 of 40 units, random selection | TWFE imputation | -0.10 (0.05) | 0.75 | 0.890 | 0.870 |
| | TROP local | -0.04 (0.02) | 0.30 | 0.940 | 0.890 |
| | TROP global | -0.03 (0.02) | 0.29 | 0.925 | 0.930 |
| A2: 30 of 80 units, random selection | TWFE imputation | -0.08 (0.03) | 0.31 | 0.920 | 0.930 |
| | TROP local | -0.04 (0.02) | 0.17 | 0.970 | 0.960 |
| | TROP global | -0.03 (0.02) | 0.16 | 0.990 | 0.980 |
| B: 10 of 40 units, selection on the trend | TWFE imputation | 1.20 (0.05) | 1.37 | 0.500 | 0.500 |
| | TROP local | 0.64 (0.02) | 0.72 | 0.475 | 0.390 |
| | TROP global | 0.37 (0.02) | 0.46 | 0.700 | 0.645 |

With random selection, the bias of each estimator is less than three Monte
Carlo standard errors. The RMSE of TROP is 44% to 62% less than the RMSE of
TWFE imputation. All three methods use the same panels. Thus, the small
negative biases are the same noise for all methods.

For TROP, the coverage of the normal interval is near 95%. With 100 draws, the
coverage of the percentile interval is only 89% for the local method in
simulation A. For this reason, the normal interval is the default. In
simulation A2 with 30 treated units, no coverage is less than 95% by more than
two Monte Carlo standard errors. The global interval has a coverage of 0.99,
so it is wider than necessary. The paper also says that the bootstrap is
correct only with many treated units.

In simulation B, the selection of treated units uses their loading on the
trend factor. Thus, TWFE imputation has a large bias. TROP removes
approximately half of this bias with the local method and two thirds with the
global method. But some bias stays, and the coverage of the intervals is less
than 95%.

This result gives no general rule about which weights can remove a bias.
Under the assumptions of the paper, a model without bias can also remove it.
Thus, this result applies to this simulation only.

## Simulations with a known effect in a staggered design

The script `staggered.R` runs the simulations C, D, and E. In a staggered
design, the treated units start the treatment in different periods. A
conventional TWFE regression can give negative weights to some effects in this
design. [Hull (2026)](https://arxiv.org/abs/2610.39995) compares this regression
with two modern estimators in real studies. The theory of the paper is for a
block design. Thus, these simulations do not test the theory.

Each panel has 60 units and 16 periods. Four cohorts of 12 units start the
treatment in periods 5, 8, 11, and 14, and stay treated. The 12 other units are
never treated. Thus, each panel has 360 treated cells. Each simulation has 200
replications. The truth is the mean of the effects in the treated cells.

- Simulation C: the effect is 1 in each treated cell. The untreated outcomes
  have unit effects, time effects, and noise. The simulation is a control for D
  and E.
- Simulation D: the effect is `delta + 0.1 * e`. The value `e` is the number of
  periods since the start of the treatment. The value `delta` is 0.8, 0.6, 0.4,
  and 0.2 for the four cohorts. The effect is between 0.2 and 1.9, and its mean
  is 1. The untreated outcomes are as in simulation C.
- Simulation E: the effect is as in D. The untreated outcomes also have the two
  factors of simulation A. The units with the largest loading on the trend
  start first, and the 12 units with the smallest loadings are never treated.
  Thus, the trends of the cohorts differ.

Three methods run on the same panels. The script stops if one of them skips a
treated cell. The global method works only for a block design, so the script
does not use it.

- `twfe_regression` is the coefficient of the treatment indicator in a
  regression with unit effects and period effects. The regression uses all 960
  cells. The standard error uses one cluster for each unit, with the factor
  G/(G-1). This method has a normal interval and no percentile interval.
- `twfe_imputation` is `lambda = c(time = 0, unit = 0, nn = Inf)`, with 200
  bootstrap draws. The fit uses only the untreated cells. The estimate is the
  mean of the effects in the 360 treated cells.
- `local` is the local method, with 100 bootstrap draws. The grid for `time` is
  0, 0.25, 0.5, and 1. The grid for `unit` is 0, 0.5, 1, and 2. The grid for
  `nn` is 0.1, 0.5, 2, and `Inf`. The leave-one-out criterion on 80 sampled
  untreated cells selects the tuning parameters. The tuning parameters stay
  fixed in the bootstrap.

The script writes three kinds of files. The files `staggered-draws-*.csv` have
the result of each method in each replication. The files
`staggered-summary-*.csv` have the summary of each simulation. The file
`staggered-twfe-weights.csv` has the weight of the regression on each group of
treated cells. A weight of 1 is the equal weight of imputation.

To run one simulation, use this command in the package folder. The arguments
are the number of replications, the number of threads, and the simulation (C,
D, or E):

```
Rscript inst/validation/staggered.R 200 8 D
```

The next table has the results. In each simulation, the three methods use the
same 200 panels.

| Simulation | Method | Bias (Monte Carlo SE) | RMSE | Coverage, normal interval | Coverage, percentile interval |
|---|---|---|---|---|---|
| C: one effect, parallel trends | TWFE regression | -0.005 (0.004) | 0.055 | 0.955 | None |
| | TWFE imputation | -0.005 (0.004) | 0.063 | 0.975 | 0.950 |
| | TROP local | -0.006 (0.005) | 0.066 | 0.970 | 0.955 |
| D: effects differ, parallel trends | TWFE regression | -0.361 (0.004) | 0.366 | 0.000 | None |
| | TWFE imputation | -0.005 (0.004) | 0.063 | 0.985 | 0.985 |
| | TROP local | -0.006 (0.005) | 0.066 | 0.980 | 0.960 |
| E: effects differ, trends differ | TWFE regression | 0.225 (0.011) | 0.277 | 0.915 | None |
| | TWFE imputation | 1.257 (0.017) | 1.280 | 0.000 | 0.005 |
| | TROP local | 1.006 (0.015) | 1.029 | 0.005 | 0.005 |

In simulation C, the three estimators have small biases. The regression has the
smallest standard deviation. It uses the treated cells in the fit, and its model
is correct when the effect is the same in each cell.

In simulation D, the regression has a large bias. Its weights give a weighted
mean of the effects of 0.64, and the true mean is 1. Three of the 30 groups of
treated cells have a negative weight. TWFE imputation and TROP have the same
biases as in simulation C. They use only the untreated cells in the fit, and
the panels have the same untreated outcomes.

In simulation D, the mean standard error of TWFE imputation is 0.080, but the
standard deviation of its estimates is 0.063. Thus, its intervals are wider
than necessary. We think that the bootstrap causes this result. It resamples
the treated units, so the number of units in each cohort changes in each draw.
The effect differs between the cohorts, so the mean effect also changes. We did
not test this cause in a saved simulation.

In simulation E, the biases of TWFE imputation and TROP are large: 1.257 and
1.006. The true effect is 1, so the mean estimate of TROP is 2.0 times the true
effect. TROP has a bias 20% smaller and an RMSE 20% smaller than TWFE
imputation. In each of the 200 panels, the estimate of TROP is closer to the
truth than the estimate of TWFE imputation. The intervals of both methods
contain the true effect in at most 1 of the 200 panels.

The cross-validation of the local method also chooses different tuning
parameters in simulation E. In simulations C and D, it chooses `nn = Inf` (no
low-rank term) in 158 of the 200 panels. In simulation E, it chooses a finite
`nn` in all 200 panels: 0.1 in 24 panels, 0.5 in 44 panels, and 2 in 132
panels. We did not study the cause of the bias that stays.

In simulation E, the regression has the smallest bias, 0.225. This result does
not show an advantage of the regression. Its bias has two parts. The weights
give -0.36, as in simulation D. The different trends add approximately 0.58. The
two parts have opposite signs, so they partly cancel, but in another design
they can have the same sign.

These simulations have limits. They use one design and one tuning grid, and the
noise is independent between the cells. Real panels can have correlation over
time. The simulations do not test the other modern estimator of Hull, which
averages the DiD estimates of each cohort and period. They also do not compare
TROP with a modern estimator on real data.

## Speed

The folder `inst/bench` has the speed benchmarks. Refer to
`inst/bench/README.md`.

## Limits of the claims

The equations of the package come from arXiv version 3 of the paper. The
documentation gives the triple robustness result only for a block design. The
package gives a warning message for a design that is not a block design. It
also gives a warning message for a bootstrap with fewer than five treated
units.

The package does not use code from the Stata ports of TROP. Some of these
ports have the GPL-3.0 license. The documentation also says that the authors
of the method did not write the package.
