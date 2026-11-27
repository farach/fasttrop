# Make the treatment and the target from start dates

`trop_rollout()` makes the treatment indicator, the target, and the
event time for a panel in long format. Each unit has one start: the
first period of the treatment. A unit without a start is never treated.
The treatment does not stop after it starts.

## Usage

``` r
trop_rollout(unit, time, start, anticipation = 0L, window = c(0, Inf))
```

## Arguments

- unit, time:

  Vectors with one value for each row of the panel. A unit must have
  only one row for each period.

- start:

  A vector with the start of each row, in the units of `time`. The value
  must be the same for all rows of a unit. Use `NA` for a unit that is
  never treated, or that starts after the last period of the panel. Each
  other value must be one of the values of `time`.

- anticipation:

  The number of periods before the start that are treated.

- window:

  Two event times. The target contains the treated cells with an event
  time from `window[1]` to `window[2]`. The default contains all treated
  periods from the start.

## Value

A data frame with one row for each row of the input and these columns:
`treated` (0 or 1), `target` (`TRUE` or `FALSE`), and `event_time` (the
number of periods since the start, or `NA` for a unit that is never
treated).

## Details

With `anticipation = k`, the `k` periods before the start are also
treated. Thus, the fit does not use them as untreated cells. The event
time is the number of periods since the true start. Thus, the
anticipation periods have the event times `-k` to `-1`.

The target contains the treated cells with an event time in `window`.
With the default window, the anticipation periods are not in the target.
Thus, they are not in the estimate. With a negative `window[1]`, the
anticipation periods in the window are in the target. Use this for
placebo estimates of the periods before the start.

Give the columns to
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) with the
arguments `treatment`, `target`, and `event_time`. With `event_time`,
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
and
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
use these event times. Then the groups are correct also when the panel
has no row for some periods, or when the anticipation periods start
before the panel. Do not give `anticipation` to these functions, because
the event times already count the anticipation periods. If you do not
give `event_time` to
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md), give
the same `anticipation` to
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
and
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md).

The periods are the sorted values of `time`. The event time counts these
periods, not the difference of the time values. Thus, a panel with one
row for each month gives the event time in months.

## See also

[`trop()`](https://farach.github.io/fasttrop/reference/trop.md),
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md),
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)

## Examples

``` r
panel <- expand.grid(unit = 1:4, month = 1:6)
starts <- c(3, 5, NA, NA)
panel$start <- starts[panel$unit]
r <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 1,
                  window = c(0, 1))
cbind(panel, r)[panel$unit == 1, ]
#>    unit month start treated target event_time
#> 1     1     1     3       0  FALSE         -2
#> 5     1     2     3       1  FALSE         -1
#> 9     1     3     3       1   TRUE          0
#> 13    1     4     3       1   TRUE          1
#> 17    1     5     3       1  FALSE          2
#> 21    1     6     3       1  FALSE          3
# Then, for example:
# trop(data, "y", "treated", "unit", "month", target = "target",
#      event_time = "event_time")
```
