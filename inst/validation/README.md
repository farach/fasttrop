# Validation

This folder has the scripts and the results for the validation of fasttrop.
The validation follows the conditions C1 to C5 of the research council that
selected the method. The package README gives a summary.

## C1. Agreement with the reference implementations

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

Each panel has five sets of tuning parameters. The sets are TWFE, weighted
TWFE, matrix completion, and two sets for the full TROP estimator. The script
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

## C2. Tuning and inference

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
  both methods. The tuning parameters stay fixed, as in Algorithm 3 of the
  paper.
- The stratified bootstrap draws treated and untreated units separately. The
  option `resample = "controls"` keeps the treated units fixed.
- With a seed, the results are the same in each run. The functions do not
  change the state of the random number generator of the caller. The number
  of threads does not change the results.

The file `tests/testthat/test-estimator.R` compares the local weights with an
independent R calculation of the formula in the paper. Another test makes sure
that `lambda = c(time = 0, unit = 0, nn = Inf)` gives the same result as
imputation with two-way fixed effects. This test has cells with `NA` outcomes.

## C3. Simulations with a known effect

The script `known_truth.R` runs the simulations A, A2, and B. The top of the
script gives the details of each simulation. The script
`summarize_known_truth.R` makes the files `known-truth-summary-*.csv` from the
saved draws. The package README gives the results and their interpretation.

## C4. Speed

The folder `inst/bench` has the speed benchmarks. Refer to
`inst/bench/README.md`.

## C5. Limits of the claims

The equations of the package come from arXiv version 3 of the paper. The
documentation gives the triple robustness result only for a block design. The
package gives a warning message for a design that is not a block design. It
also gives a warning message for a bootstrap with fewer than five treated
units.

The package does not use code from the Stata ports of TROP. Some of these
ports have the GPL-3.0 license. The documentation also says that the authors
of the method did not write the package.
