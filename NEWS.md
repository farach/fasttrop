# fasttrop 0.1.0

* First version. Implements the Triply Robust Panel (TROP) estimator of Athey,
  Imbens, Qu, and Viviano (2026): the paper's per-cell (local) estimator and the
  block (global) estimator of the authors' reference software, leave-one-out and
  placebo cross-validation, and the unit bootstrap of Algorithm 3 with normal
  (default) or percentile intervals.
* The computational core is C with OpenMP and R's BLAS/LAPACK. Results do not
  depend on the number of threads.
* On Windows, kernel threads opt out of power throttling while they run
  (`trop_control(throttle = FALSE)`), because Windows slows long-running
  background processes and this made identical native work about 1.3 times
  slower inside an R session than in a fresh process.
