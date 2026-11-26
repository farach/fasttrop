# Package index

## Estimate effects

Calculate the TROP estimate for a data frame or for matrices.

- [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) :
  Triply Robust Panel (TROP) estimator
- [`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
  : TROP estimator for outcome and treatment matrices

## Read the results

Show a fit and get its parts.

- [`print(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`summary(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`print(`*`<summary.trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`coef(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`confint(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`vcov(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  [`as.data.frame(`*`<trop>`*`)`](https://farach.github.io/fasttrop/reference/trop-methods.md)
  :

  Methods for `trop` objects

## Group effects and many tests

Get the effect of groups of cells and adjust the p-values of many tests.

- [`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
  [`print(`*`<trop_effects>`*`)`](https://farach.github.io/fasttrop/reference/trop_effects.md)
  [`as.data.frame(`*`<trop_effects>`*`)`](https://farach.github.io/fasttrop/reference/trop_effects.md)
  : Effects by event time, cohort, period, or unit
- [`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)
  : Romano-Wolf stepdown p-values
- [`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md)
  : Make bootstrap rows

## Settings

Set the solver options, the tuning grids, and the number of threads.

- [`trop_control()`](https://farach.github.io/fasttrop/reference/trop_control.md)
  : Solver options
- [`trop_grid()`](https://farach.github.io/fasttrop/reference/trop_grid.md)
  : Tuning grids
- [`trop_grid_scaled()`](https://farach.github.io/fasttrop/reference/trop_grid_scaled.md)
  : Scale the tuning grid to the data
- [`trop_threads()`](https://farach.github.io/fasttrop/reference/trop_threads.md)
  : Default number of threads

## Practice data

Make a panel with a known effect.

- [`trop_simulate()`](https://farach.github.io/fasttrop/reference/trop_simulate.md)
  : Simulated panel with a factor structure
