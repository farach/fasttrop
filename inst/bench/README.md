# Benchmarks

This folder has several benchmarks. The first benchmark measures the cost of the R
code. It compares the R package with a C program that uses the same C core.
The second benchmark compares fasttrop with the other implementations of the
TROP estimator.

## R and C with the same C core

The C program `trop_driver.c` contains `src/trop_core.c` and `src/trop_os.c`.
It runs the same sequence of C calculations as the R function `trop_run()`:
the tuning search, the estimate, and the bootstrap. The script
`run_benchmark.R` compiles the program with the compiler, flags, BLAS, and
LAPACK of R. The flags are `R CMD config CFLAGS` and `-DNDEBUG`, as in
`R CMD INSTALL`.

The script makes one plan in R and writes it to a binary file. The plan
contains the panel, the sampled cells, the placebo sets, and the bootstrap
rows. In each round, the script runs three paths on this plan in a random
sequence:

- The C program
- The R function `trop_run()`, which starts the C calculations for a prepared
  plan
- The R function `trop_matrix()`, which also examines the inputs and makes the
  random draws

The script makes sure that the C program and R give the same tuning
parameters, estimates, and bootstrap draws. It also makes sure that 1 thread
and 4 threads give the same results. The timer on Windows is
QueryPerformanceCounter, because `proc.time()` and the `omp_get_wtime()`
function of MinGW have a resolution of approximately one millisecond or worse.

To run the benchmark, use this command in the package folder:

```
Rscript inst/bench/run_benchmark.R [rounds] [threads] [inner]
```

The script writes these files:

- `results.csv`: the medians for each workload and number of threads
- `results-rounds.csv`: the times for each round
- `results-pooled.csv`: the geometric mean of the ratio of R time to C time
  for runs of one second or more, with a 95% interval

In the last run, the script used 3 rounds, 1 thread and 4 threads, and four
workloads. The times were from approximately 0.006 s to 15 s. The tuning
parameters, the estimates, and the bootstrap draws were the same for all paths
and for the two numbers of threads.

There were 18 pairs of runs of one second or more. For these pairs, the
geometric mean of the ratio was 0.996 for `trop_run()`, with a 95% interval of
0.970 to 1.023. For `trop_matrix()`, it was 1.012, with a 95% interval of
0.983 to 1.041. For the smallest workload, `trop_run()` added approximately
0.3 ms for each call. The full `trop_matrix()` call added 0.8 ms to 1.7 ms. In
the limits of this measurement, the R package has the same speed as the C
program.

| Workload | Threads | C program (s) | `trop_run()` (s) | `trop_matrix()` (s) |
|---|---|---|---|---|
| Local, fixed lambda, no bootstrap, 20 x 12 | 1 | 0.0204 | 0.0207 | 0.0221 |
| | 4 | 0.0056 | 0.0060 | 0.0064 |
| Local, fixed lambda, 100 bootstrap draws, 30 x 20 | 1 | 14.38 | 14.10 | 15.70 |
| | 4 | 3.99 | 3.88 | 3.75 |
| Local, CV on 60 cells, 50 draws, 30 x 20 | 1 | 4.80 | 4.34 | 4.69 |
| | 4 | 1.12 | 1.17 | 1.15 |
| Global, placebo CV with 50 sets, 200 draws, 50 x 30 | 1 | 7.03 | 7.21 | 7.81 |
| | 4 | 1.86 | 1.85 | 1.89 |

The values are medians of three rounds. The function `trop_run()` starts the
C calculations for a prepared plan. The function `trop_matrix()` is the full
function call, with the input checks and the random draws.

### Windows power throttling

Windows can use power throttling (EcoQoS) on a process that runs in the
background for a long time. Power throttling makes the process slower. An R
session with a long calculation can be such a process. A short C program
usually stops before Windows starts power throttling. Thus, the threads of
fasttrop tell Windows not to use power throttling while the C code runs.

When the C code stops, Windows controls the threads again. With
`trop_control(throttle = TRUE)`, Windows uses its usual policy. This option
does not change the results. On other operating systems, it has no effect.

## Comparison with the other implementations

The scripts `reference_timings.R` and `reference_timings.py` use the same two
panels and the same fixed tuning parameters. Each package uses its default
solver options. The scripts measure three tasks: the local ATT, the global ATT,
and one calculation of the leave-one-out criterion for all untreated cells.

To run the comparison, do these steps in the package folder:

1. Run `Rscript inst/bench/reference_timings.R`.
2. Clone the repositories of `trop` and `diff-diff`.
3. Run `python inst/bench/reference_timings.py <trop> <diff-diff>` with the
   folders of the two repositories.
4. Run `Rscript inst/bench/reference_timings.R combine`.

The file `reference-timings.csv` has the results. A column gives each time as
a multiple of the time of fasttrop with one thread. The panels have 20 units
and 15 periods, and 40 units and 20 periods. The tuning parameters are
(time, unit, nn) = (0.3, 0.8, 0.5). An empty cell shows that the package does
not have the method.

| Task | Panel | fasttrop, 1 thread | fasttrop, 4 threads | Authors' `trop` | `diff-diff` | `cfcompare` |
|---|---|---|---|---|---|---|
| Local ATT | 20 x 15 | 0.044 s | 0.016 s | | 3.42 s | 5.93 s |
| Local ATT | 40 x 20 | 0.33 s | 0.088 s | | 9.33 s | 49.7 s |
| Global ATT | 20 x 15 | 0.0032 s | 0.0040 s | 0.100 s | 0.551 s | |
| Global ATT | 40 x 20 | 0.010 s | 0.0084 s | 0.252 s | 0.951 s | |
| One evaluation of the leave-one-out criterion, 284 cells | 20 x 15 | 0.85 s | 0.20 s | | 56.6 s | |

With one thread, fasttrop was 25 to 31 times faster than the authors' package
for the global method. It was 29 to 174 times faster than `diff-diff` in
Python, and 136 to 152 times faster than `cfcompare`. The ATTs of fasttrop and
the Python packages had a difference of less than 1e-5. The ATTs of
`cfcompare` had a difference of up to 0.006. Thus, `cfcompare` does not
calculate exactly the same value.

`diff-diff` ran in Python only, because its optional Rust code was not
available. With the Rust code, its times can be shorter.

## The time cost of covariates

The script `covariate_timings.R` fits panels of 100 x 30 and 400 x 50 with 0,
2, and 5 covariates. It uses fixed tuning parameters, no bootstrap, and 8
threads. The file `covariate-timings.csv` has the results. The column
`ratio_to_no_covariates` gives the time with covariates divided by the time
without covariates.

The local method on the 400 x 50 panel uses 1 run only.
Some ratios are less than 1 because the solver iterations and the run times
change. The ratios show that the cost of the covariates is small compared to
the cost of the fit. To run it, use this command in the package folder:

```
Rscript inst/bench/covariate_timings.R 8 3
```

## Notes on the times

The times depend on the computer, the compiler, the BLAS library, and the
number of threads. Other work on the computer also changes the times. All
benchmarks used the reference BLAS and LAPACK of R. A faster BLAS, for example
OpenBLAS, can make the times of fasttrop shorter. When R uses a faster BLAS,
fasttrop also uses it. In the Python runs, each package used its default
solver options.

## Local time-weight cache

The script `time_cache_timings.R` compares two installed versions on the same
local cross-validation searches. The baseline is commit `e85d03f`, the squash
merge of PR #4. The new version uses a cache for time weights.

Each distinct time penalty has one vector of `2 * T - 1` values. A target
period selects a contiguous view of this vector. Thus, the cache uses space
proportional to the number of periods and distinct time penalties.
No T-by-T matrix is necessary.

The search owns the cache. All workers read it without changes.
The local estimate and bootstrap use the same representation.

Direct exponentials set the values. Repeated multiplication can retain a
positive subnormal value after the direct formula gives zero. This changes
the graph of positively weighted control cells. A regression test uses a
three-unit panel with a gap of 2,000 periods. The corrected estimator skips
the disconnected target and does not report an effect.

To compare two installed versions, run this command from the new source folder:

```
Rscript inst/bench/time_cache_timings.R BASE_LIB NEW_LIB 7
```

Each pair uses separate R processes for the two versions. Their sequence
alternates between pairs. Each worker first runs an untimed search.
The C timer excludes process startup, package loading, and panel preparation.
The short workload uses the mean of 25 searches in each round.
The other workload uses one search in each round.

The files `time-cache-rounds.csv` and `time-cache-summary.csv` contain the
measurements. The file `time-cache-session.txt` records the R environment.
Both versions used GCC 14.3.0, `-O2`, OpenMP, and the reference BLAS of R.

The grid has time and unit penalties of 0, 0.5, and 1.
The 30-by-20 panel also uses nuclear-norm penalties of 0.1, 1, and infinity.
The 12-by-120 panel disables the low-rank term.
Both searches use 60 sampled control cells and fixed seeds.

| Panel | Threads | Baseline median (s) | Cache median (s) | Paired speed ratio, 95% interval |
|---|---|---|---|---|
| 30 x 20 | 1 | 2.595 | 2.610 | 1.001, 0.965 to 1.039 |
| 30 x 20 | 4 | 0.6035 | 0.6235 | 0.975, 0.910 to 1.044 |
| 12 x 120, no low-rank term | 1 | 0.009857 | 0.009983 | 0.987, 0.919 to 1.061 |
| 12 x 120, no low-rank term | 4 | 0.005030 | 0.005239 | 0.903, 0.752 to 1.084 |

The ratios use the geometric mean of paired baseline times divided by cache
times. Values above one mean that the cache is faster.
The intervals use a Student t interval on the log ratios across seven pairs.
Every interval contains one. Thus, these measurements show no clear speed gain.

The cache median is slightly slower on the representative 30-by-20 workload.
The implementation retains the cache to correct underflow and remove repeated
time-vector calculations. It does not claim a speed increase on these panels.

The tuning parameters, cell effects, and four bootstrap draws agree between
versions. The largest CV-score difference is `1.42e-14`.
Each version gives exactly the same bits with one and four threads.
These comparisons apply to the measured panels, not all possible inputs.
