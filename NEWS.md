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

- All calculations run in C with OpenMP and with the BLAS and LAPACK libraries
  of R. The number of threads does not change the results.
- On Windows, the C code tells Windows not to use power throttling on its
  threads while it runs. In a long R session, power throttling made the same C
  code approximately 1.3 times slower. To keep the usual Windows policy, use
  `trop_control(throttle = TRUE)`.
- The function `trop_threads()` gives at most 2 threads when `R CMD check`
  sets the environment variable `_R_CHECK_LIMIT_CORES_`.
