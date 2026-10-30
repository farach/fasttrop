# Benchmarks

This folder has two benchmarks. The first benchmark measures the cost of the R
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
0.983 to 1.041. For the smallest workload, the R code added approximately
0.3 ms for each call. The package README has the full table.

### Windows power throttling

In an earlier run of this benchmark, the R path was approximately 40% slower
than the C program. For the bootstrap on one thread, R took 29.5 s and C took
21.1 s. At that time, the package did not have the change for power
throttling. A test with a small DLL found the cause. The DLL calls `dgesdd`
and `dgemm` with the same inputs in each process.

The same machine code had full speed in a new process, in a small C host
program, and in x64 Python. In a long R session, the median time was 1.34
times longer. R was slower in 28 of 30 pairs with a random sequence, and the
Wilcoxon test gave p = 2e-7.

The tests found that these items did not cause the slowdown:

- The x87 control word and the SSE control word
- The heap type of R
- The load addresses of the DLLs
- More threads that do not run

The cause was Windows power throttling (EcoQoS). When the process or the
thread was not under power throttling, the slowdown stopped. In one long R
session, the same call took 0.70 s without power throttling. With power
throttling, the median time was 0.81 s to 1.04 s.

For this reason, the threads of fasttrop tell Windows not to use power
throttling while the C code runs. When the C code stops, Windows controls the
threads again. With `trop_control(throttle = TRUE)`, Windows uses its usual
policy. This option does not change the results.

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
a multiple of the time of fasttrop with one thread.

With one thread, fasttrop was 25 to 31 times faster than the authors' package
for the global method. It was 29 to 174 times faster than `diff-diff` in
Python, and 136 to 152 times faster than `cfcompare`. The ATTs of fasttrop and
the Python packages had a difference of less than 1e-5. The ATTs of
`cfcompare` had a difference of up to 0.006. Thus, `cfcompare` does not
calculate exactly the same value.

## Test computer

The test computer was a laptop with a Snapdragon X Elite X1E80100 processor
with 12 cores and Windows 11. R was x64 R 4.6.1 under emulation on Windows on
ARM. The compiler was GCC 14.3.0 from Rtools 4.5, with the reference BLAS and
LAPACK of R. Python was x64 Python 3 with NumPy 2.5, CVXPY 1.9 with CLARABEL,
and `diff-diff` in Python only. The laptop also did other work during the
tests. Thus, all times contain some noise.
