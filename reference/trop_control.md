# Solver options

`trop_control()` sets the options of the weighted low-rank solver in C.
The estimate, the cross-validation, and the bootstrap all use this
solver.

## Usage

``` r
trop_control(
  max_iter = 100L,
  max_inner = 10L,
  tol = 1e-06,
  threads = trop_threads(),
  warm = TRUE,
  solver = c("prox", "block"),
  throttle = FALSE
)
```

## Arguments

- max_iter:

  The maximum number of outer iterations. With `solver = "block"`, each
  outer iteration is one sweep of the blocks. With `solver = "prox"`,
  the maximum number of proximal steps is `max_iter * max_inner`.

- max_inner:

  The number of proximal steps in each outer iteration. Each step is a
  step of FISTA, an accelerated proximal gradient method.

- tol:

  The convergence tolerance. The solver stops when the largest absolute
  change of the iterates between two steps is less than `tol`.

- threads:

  The number of 'OpenMP' threads. The default comes from
  [`trop_threads()`](https://farach.github.io/fasttrop/reference/trop_threads.md).

- warm:

  Logical. If `TRUE`, the fit for a cell starts from the solution for
  the previous cell of the same unit. If `FALSE`, each fit starts from
  zero. Both options give the same results in the limits of `tol`. The
  number of threads does not change the results.

- solver:

  Use `"prox"`, the default, or `"block"`. `"prox"` is an accelerated
  proximal gradient method on the low-rank term. At each step, it
  calculates the unit and time effects again, and it uses an adaptive
  restart. `"block"` is the block coordinate method of the reference
  implementations. Both methods find the minimum of the same convex
  objective. `"prox"` usually calculates fewer singular value
  decompositions.

- throttle:

  Logical. Windows can use power throttling (EcoQoS) to make processes
  that run in the background for a long time slower. An R session with a
  long calculation can be such a process.

  With `FALSE`, the default, each thread that runs the C code of
  fasttrop tells Windows not to use power throttling on it. When the C
  code stops, Windows controls the thread again. With `TRUE`, Windows
  uses its usual policy. On other operating systems, this option has no
  effect.

## Value

A list of class `trop_control`.

## Examples

``` r
trop_control(tol = 1e-8, threads = 2)
#> $max_iter
#> [1] 100
#> 
#> $max_inner
#> [1] 10
#> 
#> $tol
#> [1] 1e-08
#> 
#> $threads
#> [1] 2
#> 
#> $warm
#> [1] 1
#> 
#> $solver
#> [1] 0
#> 
#> $throttle
#> [1] 0
#> 
#> attr(,"class")
#> [1] "trop_control"
```
