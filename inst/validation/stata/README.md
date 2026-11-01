# Stata spot check for covariates

These scripts compare fasttrop with the Stata package `trop` from
gorgeousfish/TROP. The scripts only call that package. They do not contain
code from it. That package has the AGPL license.

## Steps

1. In Stata, install the package:

   ```stata
   net install trop, from("https://raw.githubusercontent.com/gorgeousfish/TROP/main") replace
   ```

2. In Stata, go to the package folder of fasttrop and run the do-file:

   ```stata
   cd "path/to/fasttrop"
   do inst/validation/stata/trop_covariates.do
   ```

   The do-file reads the test panels `tests/testthat/fixtures/XC*.csv`. It
   fits each panel with five sets of fixed tuning parameters, with the local
   method and the global method, without and with covariates. It writes
   `inst/validation/stata/stata-covariates.csv`.

3. In R, in the package folder, run:

   ```sh
   Rscript inst/validation/stata/compare_stata.R
   ```

   The script fits the same cases with fasttrop. It writes
   `inst/validation/stata/stata-comparison.csv`.

## How to read the comparison

- `att_diff` is the difference between the ATTs of the two packages.
- `shift_diff` is the difference between the changes in the ATT when you add
  the covariates. A small `shift_diff` shows that the two packages use the
  covariates in the same way, also when their solvers differ.
- The Stata package gives an error for the global method with a staggered
  design. These rows have a value of `rc` that is not 0.
- The Stata package uses `1e10` for a tuning parameter of infinity.
