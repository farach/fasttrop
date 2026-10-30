# What TROP can and cannot do

Key points

- The authors suggest TROP or synthetic difference in differences as a
  simple choice for panel data. No method works for every data set.
- The triple robustness result has a proof only for a block design.
- The bootstrap interval is reliable only with many treated units.
- In our simulation with larger trends for the treated units, TROP
  reduced the bias, but some bias stayed.
- fasttrop cannot use covariates.

## What TROP can do

- Estimate the average treatment effect on the treated (ATT) for block,
  staggered, and general designs, with the local method.
- Use panels in which some outcomes are not in the data. The fits do not
  use these cells.
- Adapt to the data. Cross-validation uses the data to set the tuning
  parameters. They control the unit weights, the time weights, and the
  low-rank term.
- Give a standard error and a confidence interval with a unit bootstrap.

## The limits of the method

### The triple robustness result has a proof only for a block design

The paper gives a limit on the bias of TROP (Theorem 5.1). This limit is
the reason for the name “triply robust”. The proof is only for a block
design, in which all treated units start the treatment in the same
period and stay treated.

fasttrop also calculates estimates for staggered and general designs,
but it gives a warning message. For these designs, the estimate is
available, but the triple robustness result has no proof. The proof also
uses fixed weights, and TROP selects its weights from the data.

### The bootstrap interval is reliable only with many treated units

The paper says that its bootstrap is correct only with many treated
units. fasttrop gives a warning message for a panel with fewer than five
treated units. Five or more treated units do not make the interval
reliable. With one or two treated units, for example one state with a
new law, the interval can be wrong.

In our simulation A, a random draw selected 10 treated units. The normal
intervals of TROP had a coverage of 92.5% to 94%, near the target of
95%. In simulation B, the selection used the trend of the units. There,
the same intervals had a coverage of only 47.5% to 70%.

fasttrop also keeps the tuning parameters fixed in the bootstrap. Thus,
the interval does not include the variation that comes from the
selection of the tuning parameters.

### Different trends can leave some bias

TROP calculates the outcomes without the treatment from the untreated
cells. In our simulation B, the selection favored units with larger
trends. TROP reduced the bias, but some bias stayed. This result gives
no general rule about which weights can remove a bias. Under the
assumptions of the paper, a model without bias can also remove it.

In simulation A, a random draw selects the treated units. In simulation
B, the selection favors the units with the largest trends. The true
effect is 1 in both simulations.

![A dot chart in two panels. In simulation A, the bias is -0.10 for
TWFE, -0.04 for TROP local, and -0.03 for TROP global. In simulation B,
the bias is 1.20 for TWFE, 0.64 for TROP local, and 0.37 for TROP
global.](limits_files/figure-html/bias-fig-1.png)

Bias of three estimators, with 200 simulated panels in each simulation.
In simulation A, a random draw selects the treated units. In simulation
B, the selection favors the units with the largest trends. The bars show
two Monte Carlo standard errors.

In simulation A, each bias is small, and each one is less than three
Monte Carlo standard errors from zero. In simulation B, TWFE has a large
bias. TROP removes approximately half of this bias with the local method
and two thirds with the global method, but some bias stays. The
intervals in simulation B also have a coverage below 95%. The article
[How we checked
fasttrop](https://farach.github.io/fasttrop/articles/validation.md)
shows the saved results for TWFE and the two TROP methods.

### General designs assume no carryover

In a general design, the treatment stops for some units. The estimator
then uses the cells after the treatment as untreated cells. This is
correct only if the treatment has no effect after it stops. The
estimator also assumes that the effect does not change with the earlier
treatment periods. fasttrop gives a warning message for a general
design.

### Assumptions in the paper

The paper assumes no spillovers and no dynamic effects. TROP uses the
cells of the treated units before the treatment as untreated cells.
Thus, the treatment must not change outcomes before it starts. The
treatment of one unit also must not change the outcomes of other units.
Outcomes can still change over time for other reasons.

## The limits of this package

These limits apply to fasttrop 0.1.0:

- The package cannot use covariates.
- The global method works only for a block design, and it must have at
  least two units that are never treated.
- Outcomes can be `NA`, but the tests for the tuning and for the
  bootstrap use only a small number of `NA` cells.
- The equations come from arXiv version 3 of the paper. A line-by-line
  comparison with the journal version is not complete.
- We ran the validation scripts and the benchmark scripts only on
  Windows. GitHub Actions runs `R CMD check` on Linux, macOS, and
  Windows.

## Good practice

1.  Make simulations that copy your data, and compare the methods on
    them. This is the preferred approach of the authors of TROP. The
    function
    [`trop_simulate()`](https://farach.github.io/fasttrop/reference/trop_simulate.md)
    makes a simple panel. For a real check, the simulation must copy the
    patterns of your own data.
2.  Compare TROP with TWFE. The article [How to read your
    results](https://farach.github.io/fasttrop/articles/read-your-results.md)
    shows how.
3.  Make sure that the design that fasttrop reports is the design that
    you expect.
4.  Use enough bootstrap draws for the final result. The default is 200.
5.  Report the method, the tuning parameters, the number of treated
    units, the number of bootstrap draws, and the type of interval.
6.  With few treated units, use the interval with care, and say so in
    your report.
