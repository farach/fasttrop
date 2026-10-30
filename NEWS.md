# fasttrop 0.1.0

This is the first version of the package.

## Estimator

- The package calculates the Triply Robust Panel (TROP) estimator of Athey,
  Imbens, Qu, and Viviano (2026).
- The local method (`method = "local"`) is the estimator of the paper for each
  treated cell. The global method (`method = "global"`) is the block estimator
  of the Python package from the authors.
- Leave-one-out cross-validation or placebo cross-validation selects the
  tuning parameters.
- The unit bootstrap of Algorithm 3 gives standard errors. The default
  confidence interval is a normal interval. A percentile interval is also
  available.

## Speed

- The model fits run in C with OpenMP and with the BLAS and LAPACK libraries
  of R. The number of threads does not change the results.
- On Windows, the C code tells Windows not to use power throttling on its
  threads while it runs. Power throttling can make a long calculation in an R
  session slower. To keep the usual Windows policy, use
  `trop_control(throttle = TRUE)`.
- The function `trop_threads()` gives at most 2 threads when `R CMD check`
  sets the environment variable `_R_CHECK_LIMIT_CORES_`.

## Documentation

- The package has a website with articles in plain words:
  <https://farach.github.io/fasttrop/>.
- The vignette "Get started with fasttrop" shows a first analysis step by
  step.
- `citation("fasttrop")` gives the references for the paper and the package.
