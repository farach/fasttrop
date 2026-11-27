# Plan, run, and examine a large study

Key points

- Use
  [`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)
  to make the treatment, the target, and the event times from the start
  date of each unit. Give the event times to
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) with
  `event_time`.
- With `event_time`, the groups are correct also when the panel has no
  row for some periods.
- Use
  [`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
  to count the cells of each group before you calculate the bootstrap.
  It does not show the effects.
- Give all sets of groups to one call of
  [`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md).
  Then the C code calculates the estimator only one time for each draw.
- Use `checkpoint` to save the draws to a file during a long run. If R
  stops, the same call continues from the file.
- Read the column `n_draws` and the diagnostics of the draws before you
  report a result.

## A practice rollout

In this practice panel, 60 units have 14 monthly periods. Twenty-four
units start the treatment in month 7, 9, or 11. The other units do not
start in the panel. The units know their start month. Thus, the outcome
changes by 0.5 one month before the start. The true effect is 1 at the
start, and it increases by 0.2 each month after the start.

``` r

library(fasttrop)
set.seed(42)
units <- sprintf("org%02d", 1:60)
start_month <- c(rep(c(7, 9, 11), each = 8), rep(NA, 36))
panel <- expand.grid(unit = units, month = 1:14, stringsAsFactors = FALSE)
panel$start <- start_month[match(panel$unit, units)]
k <- panel$month - panel$start
panel$y <- rnorm(60)[match(panel$unit, units)] + cumsum(rnorm(14, 0, 0.3))[panel$month] +
  rnorm(nrow(panel), 0, 0.3) +
  ifelse(!is.na(k) & k >= 0, 1 + 0.2 * k, 0) + ifelse(!is.na(k) & k == -1, 0.5, 0)
```

## Make the treatment, the target, and the event times

The month before the start is not a good control, because the outcome
changes in that month.
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)
marks this month as treated with `anticipation = 1`. A treated cell is
not a control for any fit.

The target is the set of treated cells in the estimate. This study
examines the first 4 months of the treatment. Thus, the window is the
event times 0 to 3. The anticipation month is treated, but it is not in
the target.

The event time is the number of months since the true start. Thus, the
anticipation month has the event time -1.

``` r

r <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 1,
                  window = c(0, 3))
panel <- cbind(panel, r)
panel[panel$unit == "org01", c("month", "start", "treated", "target", "event_time")]
#>     month start treated target event_time
#> 1       1     7       0  FALSE         -6
#> 61      2     7       0  FALSE         -5
#> 121     3     7       0  FALSE         -4
#> 181     4     7       0  FALSE         -3
#> 241     5     7       0  FALSE         -2
#> 301     6     7       1  FALSE         -1
#> 361     7     7       1   TRUE          0
#> 421     8     7       1   TRUE          1
#> 481     9     7       1   TRUE          2
#> 541    10     7       1   TRUE          3
#> 601    11     7       1  FALSE          4
#> 661    12     7       1  FALSE          5
#> 721    13     7       1  FALSE          6
#> 781    14     7       1  FALSE          7
```

## Fit one time without the bootstrap

The fit uses fixed tuning parameters to keep the build of this page
fast. In a real study, select them with cross-validation before you
examine the treated outcomes. Use `se = "none"`, because the bootstrap
comes later with all the groups.

The argument `event_time` keeps the event times of
[`trop_rollout()`](https://farach.github.io/fasttrop/reference/trop_rollout.md)
in the fit. Then
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
and
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
use them for the groups.

``` r

lambda <- c(time = 0.2, unit = 0.5, nn = 1)
fit <- trop(panel, outcome = "y", treatment = "treated", unit = "unit",
            time = "month", target = "target", event_time = "event_time",
            lambda = lambda, se = "none")
fit
#> Triply Robust Panel (TROP) estimate
#>   Method: local | Design: staggered | Units: 60 | Periods: 14 | Treated cells: 96
#>   Cells outside of `target`: 72 (not in the estimate)
#>   ATT: 1.242
#>   lambda: time = 0.2, unit = 0.5, nn = 1
```

## Count the cells before you see the results

[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
shows the cells of each group, but not the effects. Use it to make sure
that the estimate contains the correct cells for your question.

``` r

trop_cells(fit)
#> TROP cells by event time
#>   The event times come from `event_time` of the fit.
#>   Weights: each cell has the same weight.
#>  event_time n_units n_cells n_excluded n_skipped share
#>          -1       0       0         24         0  0.00
#>           0      24      24          0         0  0.25
#>           1      24      24          0         0  0.25
#>           2      24      24          0         0  0.25
#>           3      24      24          0         0  0.25
#>           4       0       0         16         0  0.00
#>           5       0       0         16         0  0.00
#>           6       0       0          8         0  0.00
#>           7       0       0          8         0  0.00
#> Cells in the estimate: 96. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 72 outside of the target.
```

Each row is an event time. The column `n_cells` counts the cells in the
estimate. The column `n_excluded` counts the treated cells outside of
the target. Event time -1 is the anticipation month. Its cells are
treated, but they are not in the estimate. The column `share` gives the
weight of each event time in the overall estimate.

## Why the fit must keep the event times

A fit without `event_time` finds the start of each run from the
treatment. In this panel, the first treated month is the anticipation
month. Thus, all labels move by one month:

``` r

fit_w <- trop(panel, outcome = "y", treatment = "treated", unit = "unit",
              time = "month", target = "target", lambda = lambda, se = "none")
head(trop_cells(fit_w), 3)
#> TROP cells by event time
#>   Weights: each cell has the same weight.
#>  event_time n_units n_cells n_excluded n_skipped share
#>           0       0       0         24         0  0.00
#>           1      24      24          0         0  0.25
#>           2      24      24          0         0  0.25
#> Cells in the estimate: 96. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 72 outside of the target.
```

With `anticipation = 1`,
[`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
and
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
move the labels back by one month. This method works only if the panel
has a row for the first treated month of each unit. In the next panel,
unit `org01` has no row for month 6. Month 6 is its anticipation month
and its first treated month. Thus, its first treated row is its true
start.

``` r

gap <- panel[!(panel$unit == "org01" & panel$month == 6), ]
fit_gap <- trop(gap, outcome = "y", treatment = "treated", unit = "unit",
                time = "month", target = "target", lambda = lambda, se = "none")
#> Warning: For 1 unit(s), the data have no rows for the periods just before the
#> first treated row of a run. Thus, the start of these runs is not certain, and
#> the function uses the first treated row. For correct event times, give
#> `event_time`, or add rows with an NA outcome. Units: org01.
#> Warning: Theorem 5.1 of the paper (triple robustness) has a proof only for a
#> block design. This design is staggered.
```

The first message tells you that the start of the treatment of `org01`
is not certain. The second message is about the design, and it is also
true for the first fit. With `anticipation = 1`, event time -1 now has a
cell in the estimate. This cell is the true start of `org01`:

``` r

head(trop_cells(fit_gap, anticipation = 1), 3)
#> TROP cells by event time
#>   Event time 0 is the first period after 1 anticipation period(s).
#>   Weights: each cell has the same weight.
#>  event_time n_units n_cells n_excluded n_skipped  share
#>          -1       1       1         23         0 0.0104
#>           0      24      24          0         0 0.2500
#>           1      24      24          0         0 0.2500
#> Cells in the estimate: 96. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 71 outside of the target.
```

With `event_time`, the groups are correct. A month without a row does
not change the event times of the other months.

``` r

fit_gap <- trop(gap, outcome = "y", treatment = "treated", unit = "unit",
                time = "month", target = "target", event_time = "event_time",
                lambda = lambda, se = "none")
head(trop_cells(fit_gap), 3)
#> TROP cells by event time
#>   The event times come from `event_time` of the fit.
#>   Weights: each cell has the same weight.
#>  event_time n_units n_cells n_excluded n_skipped share
#>          -1       0       0         23         0  0.00
#>           0      24      24          0         0  0.25
#>           1      24      24          0         0  0.25
#> Cells in the estimate: 96. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 71 outside of the target.
```

A month without a row in the middle of a treatment run does not stop the
run. [`trop()`](https://farach.github.io/fasttrop/reference/trop.md)
gives such a cell the treatment of the rows before and after it. The
cell has no outcome. Thus, it is not in the estimate, and the labels of
the other months do not change, also with `anticipation`.

## Choose the weights

By default, each cell has the same weight. Thus, a unit with more cells
in a group has more weight in that group. With `weights = "unit"`, each
unit has the same total weight in each group. In this study, each
treated unit has 4 cells in the window. Thus, the two rules give the
same estimates here.

With `unit_weights`, each unit gets its own weight, for example its
number of workers before the treatment. Use a value that the treatment
cannot change. A unit with a weight of 0 is not in the estimate. The
weights are relative. Thus, only the ratios of the weights change the
estimate.

``` r

workers <- setNames(round(runif(60, 50, 500)), units)
trop_cells(fit, by = "overall", weights = "unit", unit_weights = workers)
#> TROP cells of the overall estimate
#>   Weights: the total weight of each unit in each group is its value in `unit_weights`.
#>  overall n_units n_cells n_excluded n_skipped share
#>  overall      24      96         72         0     1
#> Cells in the estimate: 96. The column share is the weight of each group in the overall estimate.
#> Treated cells that are not in the estimate: 72 outside of the target.
```

## Run the bootstrap one time for all groups

Make the bootstrap rows one time with
[`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md).
Then give all sets of groups to one call of
[`trop_effects()`](https://farach.github.io/fasttrop/reference/trop_effects.md)
as a named list. Each draw then calculates the estimator only one time,
for all the groups.

A long bootstrap can take hours. With `checkpoint`, the function saves
the draws to a file after each block of draws. With `progress = TRUE`,
it shows the number of draws that are done and an estimate of the time
that is left. This page uses 100 draws. A real study must use more
draws.

``` r

W <- tapply(panel$treated, list(panel$unit, panel$month), sum)
rows <- trop_boot_rows(W, n_boot = 100, seed = 7)
ck <- tempfile("my-study-draws-", fileext = ".rds")
res <- trop_effects(fit, by = list(att = "overall", event = "event_time",
                                   cohort = "cohort"),
                    boot_rows = rows, checkpoint = ck, progress = TRUE)
#> Bootstrap: 100 draws with 4 thread(s). 0 draw(s) are done.
#> Bootstrap draws: 16 of 100 done. Time: 1 s. Time left: about 3 s.
#> Bootstrap draws: 32 of 100 done. Time: 1 s. Time left: about 2 s.
#> Bootstrap draws: 48 of 100 done. Time: 2 s. Time left: about 2 s.
#> Bootstrap draws: 64 of 100 done. Time: 2 s. Time left: about 1 s.
#> Bootstrap draws: 80 of 100 done. Time: 3 s. Time left: about 1 s.
#> Bootstrap draws: 96 of 100 done. Time: 3 s. Time left: about 0 s.
#> Bootstrap draws: 100 of 100 done. Time: 3 s. Time left: about 0 s.
res$event
#> TROP effects by event time
#>   The event times come from `event_time` of the fit.
#>   Weights: each cell has the same weight.
#>  event_time estimate n_cells n_units std.error conf.low conf.high band.low
#>           0   0.8926      24      24   0.07136   0.7527     1.032   0.7273
#>           1   1.1279      24      24   0.07412   0.9826     1.273   0.9563
#>           2   1.4224      24      24   0.07171   1.2819     1.563   1.2564
#>           3   1.5266      24      24   0.06964   1.3901     1.663   1.3653
#>  band.high n_draws
#>      1.058     100
#>      1.299     100
#>      1.588     100
#>      1.688     100
#> Bootstrap: 100 draws. The band is simultaneous (sup-t) at the 95% level. It uses the 100 draws with a value for all rows.
```

If R stops before the end, run the same call again. The function reads
the file and calculates only the draws that are not in it. The results
are equal to the results of one run without a stop. The file has a
fingerprint of the calculation. If you change the data, the groups, the
weights, the rows, or the tuning parameters, the function stops with an
error. It does not use the old draws.

## Read the diagnostics

The column `n_draws` gives the number of draws with a value for each
group. A group has no value in a draw if the draw has no treated cell of
the group. The band uses only the draws with a value for all groups. The
last lines of the print show this number.

Many draws can have no value for a group, for example for a small
cohort. Then you can make the draws with the cohort of each unit as
`strata` in
[`trop_boot_rows()`](https://farach.github.io/fasttrop/reference/trop_boot_rows.md).
The units that are never treated must have their own stratum, because
`strata` cannot have `NA` values. Each draw then has the same number of
units of each cohort as the data:

``` r

first <- tapply(panel$start, panel$unit, function(s) s[1])
cohort <- ifelse(is.na(first), "never", first)
rows_c <- trop_boot_rows(W, n_boot = 100, strata = cohort, seed = 7)
table(cohort[rows_c[1, ]])
#> 
#>    11     7     9 never 
#>     8     8     8    36
```

A group can still have no value in a draw if the untreated units of the
draw do not identify its cells. Thus, examine `n_draws` and the column
`n_dropped` of the diagnostics also with these draws. A cohort with few
units has little variation in such draws.

The attribute `diagnostics` counts the fits of each draw. A fit that
does not converge stops at the iteration limit. A fit that fails has no
result, because the singular value decomposition failed or the result is
not finite.

``` r

head(attr(res$event, "diagnostics"), 3)
#>   draw n_used n_skipped n_nonconverged n_failed n_dropped
#> 1    1     96         0              0        0         0
#> 2    2     96         0              0        0         0
#> 3    3     96         0              0        0         0
```

The next fit uses a very low iteration limit to show the messages. Do
not use such a limit in a real study.

``` r

low <- fit
low$control <- trop_control(max_iter = 1, max_inner = 1)
trop_effects(low, by = "overall", boot_rows = rows[1:10, ])
#> Warning: In 10 of 10 bootstrap draws, 960 fit(s) did not converge. Use a larger
#> `max_iter` in trop_control().
#> TROP effects of the overall estimate
#>   Weights: each cell has the same weight.
#>  overall estimate n_cells n_units std.error conf.low conf.high band.low
#>  overall    1.242      96      24   0.03602    1.172     1.313       NA
#>  band.high n_draws
#>         NA      10
#> Bootstrap: 10 draws. The band needs 2 or more rows.
#> Bootstrap fits: 960 fit(s) in 10 draw(s) did not converge.
#> Refer to attr(x, "diagnostics") for the counts of each draw.
```

If a draw has fits that did not converge, use a larger `max_iter` in
[`trop_control()`](https://farach.github.io/fasttrop/reference/trop_control.md).
If a draw cannot identify some cells, the draw uses fewer cells than the
estimate. Then the standard error is not the standard error of the same
estimand. Find the reason before you report the result.

## Combine the groups

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) give the estimates and
their covariance from the draws.
[`trop_contrast()`](https://farach.github.io/fasttrop/reference/trop_contrast.md)
calculates linear combinations of the estimates from the same draws. It
does not calculate the estimator again.

``` r

ct <- trop_contrast(res$event, list(
  mean_0_3 = c("0" = 1, "1" = 1, "2" = 1, "3" = 1) / 4,
  change_0_to_3 = c("3" = 1, "0" = -1)
))
ct
#> TROP effects by contrast
#>   Each row is a linear combination of estimates.
#>       contrast estimate std.error conf.low conf.high band.low band.high n_draws
#>       mean_0_3    1.242   0.04214   1.1598    1.3250   1.1544    1.3304     100
#>  change_0_to_3    0.634   0.08452   0.4684    0.7997   0.4575    0.8105     100
#> Bootstrap: 100 draws. The band is simultaneous (sup-t) at the 95% level. It uses the 100 draws with a value for all rows.
```

The result is a `trop_effects` object. Thus,
[`trop_stepdown()`](https://farach.github.io/fasttrop/reference/trop_stepdown.md)
can adjust the p-values of all the tests in a family together.

``` r

trop_stepdown(list(event = res$event, contrast = ct))
#> Romano-Wolf stepdown p-values
#>                     name estimate std.error      t    p_raw    p_adj
#>                 event: 0   0.8926   0.07136 12.507 0.009901 0.009901
#>                 event: 1   1.1279   0.07412 15.218 0.009901 0.009901
#>                 event: 2   1.4224   0.07171 19.836 0.009901 0.009901
#>                 event: 3   1.5266   0.06964 21.922 0.009901 0.009901
#>       contrast: mean_0_3   1.2424   0.04214 29.480 0.009901 0.009901
#>  contrast: change_0_to_3   0.6340   0.08452  7.501 0.009901 0.009901
#> Draws: 100 of 100 have a value for all estimates. The p-values use these draws.
```

## Placebo periods for a sensitivity analysis

The ‘HonestDiD’ package gives confidence sets when parallel trends is
not exactly true. It uses event-time estimates before and after the
start, relative to a reference period. To get estimates before the
start, mark these periods as treated and keep them in the target. Then
the fit estimates them as placebo effects.

In the next fit, the 3 months before the start are placebo periods. The
window starts at event time -3. The event times do not change, because
they count the months from the true start.

``` r

rp <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 3,
                   window = c(-3, 3))
panel$placebo <- rp$treated
panel$placebo_target <- rp$target
fit_p <- trop(panel, outcome = "y", treatment = "placebo", unit = "unit",
              time = "month", target = "placebo_target", event_time = "event_time",
              lambda = lambda, boot_rows = rows)
ev_p <- trop_effects(fit_p, by = "event_time")
ev_p[, c("event_time", "estimate", "std.error", "n_draws")]
#> TROP effects by 
#>  event_time  estimate std.error n_draws
#>          -3  0.027924   0.07457     100
#>          -2 -0.004656   0.05621     100
#>          -1  0.510645   0.07404     100
#>           0  0.895630   0.07714     100
#>           1  1.137528   0.07893     100
#>           2  1.428386   0.06883     100
#>           3  1.536789   0.07369     100
```

The estimate at event time -1 is not near 0, because the units change
one month before the start. In this design, the placebo months are not
controls. Thus, the effects after the start are different from the
effects of the first fit.

[`trop_honest_inputs()`](https://farach.github.io/fasttrop/reference/trop_honest_inputs.md)
subtracts the estimate of the reference period from the other estimates.
It also gives their covariance from the draws. The element `l_vec` gives
the weight of each event time after the reference in the result of
‘HonestDiD’.

``` r

inp <- trop_honest_inputs(ev_p, reference = -2)
round(inp$betahat, 3)
#>    -3    -1     0     1     2     3 
#> 0.033 0.515 0.900 1.142 1.433 1.541
c(pre = inp$numPrePeriods, post = inp$numPostPeriods)
#>  pre post 
#>    1    5
inp$l_vec
#> -1  0  1  2  3 
#>  0  1  0  0  0
```

With the reference -2, event time -1 is a period after the reference.
Thus, the sensitivity analysis treats the anticipation month as a
treated period. By default, ‘HonestDiD’ examines the first period after
the reference. Here, this period is event time -1. Give `l_vec` to
examine event time 0. With `event_times = 0:3`, `l_vec` gives the mean
of the event times 0 to 3.

The next code gives the inputs to ‘HonestDiD’. This website does not run
this code, because it does not install ‘HonestDiD’.

``` r

HonestDiD::createSensitivityResults_relativeMagnitudes(
  betahat = inp$betahat, sigma = inp$sigma,
  numPrePeriods = inp$numPrePeriods, numPostPeriods = inp$numPostPeriods,
  l_vec = inp$l_vec, Mbarvec = c(0.5, 1, 2)
)
```

Limits

- [`trop_cells()`](https://farach.github.io/fasttrop/reference/trop_cells.md)
  counts cells. It cannot tell you if the estimand answers your
  question. Write the estimand before you see the results.
- Without `event_time`, the groups use the first treated row of each
  run. If the panel has no row for the first treated period of a run,
  the event times and the cohort of that run are wrong.
  [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) gives
  a warning when the start of a run is not certain.
- The bootstrap keeps the tuning parameters fixed.
- A checkpoint file contains the draws, but not the panel. Keep it in a
  folder on your computer. Do not change the calculation between two
  runs with the same file.
- [`trop_honest_inputs()`](https://farach.github.io/fasttrop/reference/trop_honest_inputs.md)
  does not accept event times with gaps. ‘HonestDiD’ treats the
  estimates as consecutive periods.
- The placebo design changes the fit, because the placebo periods are
  not controls.
