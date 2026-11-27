# Inputs for a sensitivity analysis of pre-trends

`trop_honest_inputs()` makes the inputs of the 'HonestDiD' package from
event-time effects. The package gives confidence sets for the effects
when parallel trends is not exactly true.

## Usage

``` r
trop_honest_inputs(x, reference = -1L, event_times = 0L)
```

## Arguments

- x:

  A `trop_effects` object with `by = "event_time"` and a bootstrap.

- reference:

  The event time of the reference period. It must be less than 0.

- event_times:

  The event times of the parameter of interest. They must be after the
  reference. The element `l_vec` of the result gives their mean.

## Value

A list with `betahat`, `sigma`, `numPrePeriods`, and `numPostPeriods`,
the arguments of the functions of 'HonestDiD'. It also has `l_vec` (the
weights of the post periods for `event_times`), `event_time` (the event
time of each element of `betahat`), `reference`, and `n_draws` (the
number of draws in `sigma`).

## Details

'HonestDiD' uses estimates that are relative to a reference period. The
function subtracts the estimate of the `reference` event time from the
other estimates. It also calculates their covariance from the bootstrap
draws. The covariance uses only the draws with a value for all event
times.

The negative event times must be placebo estimates. To get them, set the
periods before the start to `W = 1` and keep them in the target. The
function
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)
makes such a design with a negative `window[1]`. Give its event times to
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) with
`event_time`. Without `event_time`, give the number of these periods to
`anticipation` of
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md).

In this design, the fit does not use these periods as untreated cells.
Thus, the effects after the start can be different from the effects of a
fit that uses these periods.

'HonestDiD' treats the estimates as consecutive periods. Thus, the event
times must have no gaps. The event times before the reference are the
periods before the treatment. The event times after the reference are
the periods after the treatment, with the anticipation periods after the
reference.

By default, 'HonestDiD' gives a confidence set for the first period
after the reference. With `reference = -2`, this period is the event
time -1, not the event time 0. Thus, give `l_vec = inp$l_vec` to the
functions of 'HonestDiD'. The vector `l_vec` selects the mean of the
event times in `event_times`.

## See also

[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)

## Examples

``` r
sim <- trop_simulate(n_units = 24, n_periods = 10, n_treated = 6,
                     n_post = 3, effect = 1, seed = 1)
# The start of each treated unit is its first treated period.
start <- ave(ifelse(sim$d == 1, sim$time, NA), sim$unit,
             FUN = function(v) if (all(is.na(v))) NA else min(v, na.rm = TRUE))
# Two placebo periods before the start, in the target.
sim <- cbind(sim, trop_rollout(sim$unit, sim$time, start, anticipation = 2,
                               window = c(-2, Inf)))
fit <- trop(sim, "y", "treated", "unit", "time", target = "target",
            event_time = "event_time", lambda = c(time = 0.2, unit = 0.5, nn = 1),
            se = "none", control = trop_control(threads = 1))
ev <- trop_effects(fit, by = "event_time", n_boot = 30, seed = 1)
inp <- trop_honest_inputs(ev, reference = -1)
inp$betahat
#>          -2           0           1           2 
#> -0.05585504  1.18442129  1.20402776  1.36710159 
inp$l_vec
#> 0 1 2 
#> 1 0 0 
```
