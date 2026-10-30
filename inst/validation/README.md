# Validation

This folder documents how fasttrop was checked against the conditions set by
the research council that selected the method (see the package README).

## C1. Agreement with the reference implementations

Reference implementations, pinned to these commits:

- `trop` by the paper's authors (zhaonanq/TROP-Estimator, MIT),
  commit `43738cdaa8eb8356566903bb010bbf6ce4cb31b9`. Global block estimator,
  solved with CVXPY and CLARABEL.
- `diff-diff` (igerber/diff-diff, MIT), commit
  `d46a5fd8b377923d944dc47bcf6d82233d0ebcfd`. Local (paper Eq. 12,
  Algorithm 2) and global estimators, pure Python, tolerance 1e-9.

The fixtures are four simulated panels in `tests/testthat/fixtures/`: D1
(12 units x 10 periods, block design), D2 (20 x 15, block), D3 (15 x 12,
staggered), and D4 (16 x 12, block with one treated unit). Each is fitted with
five tuning configurations: TWFE, weighted TWFE, matrix completion, and two full
TROP settings. `data-raw/oracle_fixtures.py` computes the reference answers.
fasttrop runs with tolerance 1e-10.

Largest absolute differences (`oracle_report.R`, results in
`oracle-att-comparison.csv` and `oracle-objective-comparison.csv`):

| Comparison | Largest difference |
|---|---|
| Local ATT vs diff-diff | 4.0e-8 |
| Local cell effects vs diff-diff | 1.1e-7 |
| Global ATT vs diff-diff | 2.1e-8 |
| Global cell effects vs diff-diff | 6.0e-8 |
| Global ATT vs the authors' `trop` | 2.5e-5 |
| Objective on identical weights, relative to diff-diff | 7e-15 |
| Fitted values on identical weights vs diff-diff | 2.5e-8 |

The gap to the authors' package equals the gap between the two references
themselves (their conic solver is less accurate), and is at most 3e-12 when
there is no low-rank term. On identical weights, fasttrop's objective equals
diff-diff's to machine precision, so both solvers reach the same optimum.
`tests/testthat/test-oracle.R` enforces these comparisons with tolerances
1e-6 (diff-diff), 1e-4 (authors' package), and 1e-7 relative (objective).

## C2. Tuning and inference

`tests/testthat/test-tuning-inference.R` checks the following properties.

- The leave-one-out criterion equals a direct computation in which each
  held-out cell is marked treated, so the cell is excluded from its own fit and
  from its own unit distances.
- The placebo criterion is the RMSE of placebo ATEs.
- The coordinate searches report the score of the chosen point, and the global
  search ends at a fixed point where no single grid move improves the score.
- Every bootstrap draw equals a fresh estimate on the resampled panel, for both
  methods, with tuning held fixed (paper Algorithm 3).
- Stratified resampling draws treated and untreated units separately, and the
  controls-only option keeps the treated units fixed.
- Results are reproducible with a seed, leave the caller's random state
  unchanged, and do not depend on the number of threads.

`tests/testthat/test-estimator.R` also checks the local weights against an
independent R implementation of the paper's definition, and checks that
lambda = (0, 0, Inf) reproduces two-way fixed effects imputation, including
with missing outcomes.

## C3. Known-truth Monte Carlo

`known_truth.R` runs designs A, A2, and B (described at the top of the script),
and `summarize_known_truth.R` rebuilds `known-truth-summary-*.csv` from the
saved draws. Results and their interpretation are in the package README.
