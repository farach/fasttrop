# Benchmarks

Two questions are measured here. Does calling the native core from R cost
anything compared with calling it from C? And how does fasttrop compare with
the other TROP implementations on the same computations?

## R versus C on the same core

`trop_driver.c` includes `src/trop_core.c` and `src/trop_os.c` and runs the
same sequence of native calls as `trop_run()` in R: tuning search, point
estimate, and bootstrap. `run_benchmark.R` builds it with R's own compiler,
flags (`R CMD config CFLAGS` plus `-DNDEBUG`, as `R CMD INSTALL` uses), OpenMP
flags, BLAS, and LAPACK. It then draws one plan in R (panel, sampled CV cells,
placebo sets, bootstrap rows), writes it to a binary file, and times three
paths on that plan in random order within each round:

- the C driver,
- `trop_run()`, the R function that makes the native calls for a prepared plan,
- `trop_matrix()`, the full user-facing call including validation and sampling.

The script checks that the driver and R return the same tuning parameters,
ATT, and bootstrap draws, and that results are identical with 1 and 4 threads.
Times use a high-resolution clock (QueryPerformanceCounter on Windows), because
`proc.time()` and MinGW's `omp_get_wtime()` resolve only to about a millisecond
or worse on Windows.

```
Rscript inst/bench/run_benchmark.R [rounds] [threads] [inner]
```

Results: `results.csv` (medians per configuration), `results-rounds.csv`
(every round), and `results-pooled.csv` (geometric mean of the per-round ratio
R / C for runs of at least one second, with a 95% interval).

Final run (3 rounds, 1 and 4 threads, four workloads from 0.006 s to 15 s):
across the 18 paired runs of at least one second, the geometric mean of the
time ratio R / C was 0.996 (95% interval 0.970 to 1.023) for `trop_run()` and
1.012 (0.983 to 1.041) for `trop_matrix()`. On the smallest workload the R
layer added about 0.3 ms per call. Tuning parameters, estimates, and bootstrap
draws were identical between the paths and between thread counts. The package
README has the full table.

### Windows power throttling

An earlier run of this benchmark, before the throttling change, found the R
path about 40% slower than the C driver on the same 1-thread bootstrap workload
(29.5 s against 21.1 s). A probe (identical `dgesdd` and `dgemm` calls from a
small DLL) separated the causes. The same machine code ran at full speed in a
fresh executable, in a minimal C host, and in x64 Python, but slowed by a
median factor of 1.34 inside a long-running R session (R slower in 28 of 30
randomized pairs, Wilcoxon p = 2e-7). The x87 and SSE control words were the
same in both processes, and neither R's heap type, the DLL load addresses, nor
extra idle threads reproduced the slowdown elsewhere. Opting the process or the
thread out of Windows execution-speed power throttling (EcoQoS) removed it.
Inside one long R session, the same call took 0.70 s with the opt-out and
drifted to a median of 0.81 to 1.04 s without it.

fasttrop therefore opts its kernel threads out while they run and hands the
decision back to the system afterwards. `trop_control(throttle = TRUE)` keeps
the system policy, and the setting has no effect on results.

## Matched timings against the references

`reference_timings.R` and `reference_timings.py` time the same fixed-lambda
computations on the same two panels with each package's default solver
settings: the local ATT, the global ATT, and one leave-one-out criterion
evaluation over every untreated cell.

Results are in `reference-timings.csv`, with a column giving each time as a
multiple of fasttrop's single-thread time. On one thread fasttrop was 25 to 31
times faster than the authors' package (global estimator), 29 to 174 times
faster than pure-Python diff-diff, and 136 to 152 times faster than cfcompare.
Estimates agreed with the Python packages to within 1e-5. cfcompare's differed
by up to 0.006, so it is close to but not the identical computation.

## Machine

Snapdragon X Elite X1E80100 (12 cores), Windows 11, x64 R 4.6.1 under Windows
on ARM emulation, Rtools 4.5 (GCC 14.3.0), R's reference BLAS and LAPACK, and
x64 Python 3 with NumPy 2.5, CVXPY 1.9 (CLARABEL), and pure-Python diff-diff.
The laptop was in normal use during the runs, which adds noise to every
timing.
