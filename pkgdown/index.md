# fasttrop

::: {.hero-lead}
Estimate the effect of a policy or other treatment in panel data with the
Triply Robust Panel (TROP) estimator. The model fits run in C, so the results
come fast.
:::

::: {.start-cards}
::: {.start-card}
[New here?]{.kicker}
[TROP in plain words](articles/the-big-idea.html){.stretched-link}
[What the method does and why its authors made it.]{.card-text}
:::

::: {.start-card}
[Ready to try it?]{.kicker}
[Get started](articles/fasttrop.html){.stretched-link}
[Install the package and make your first estimate.]{.card-text}
:::

::: {.start-card}
[Want the evidence?]{.kicker}
[How we checked fasttrop](articles/validation.html){.stretched-link}
[Agreement with other packages, tests, and simulations.]{.card-text}
:::
:::

## What fasttrop does

Some units get a treatment, for example a new law, and other units do not.
fasttrop estimates how the treatment changed an outcome of the treated units.
For each treated unit and period, it calculates the outcome without the
treatment from the untreated data. Three ideas control this calculation:

- Units with a similar history get more weight.
- Periods near the treated period get more weight.
- A model finds patterns that all units share over time.

The data select how much each idea counts. fasttrop then averages the effects
of all treated cells. This average is the average treatment effect on the
treated (ATT). The standard error and the confidence interval show how
uncertain the estimate is.

::: {.fact-list}
**Fast**
In our benchmarks with one thread and fixed tuning parameters, fasttrop was 25
to 174 times faster than the Python packages. Times depend on the computer.

**Checked**
On four test panels, the local ATTs of fasttrop and `diff-diff` differed by at
most 4e-8. The gap to the authors' package is up to 2.5e-5, the same as the
gap between the two Python packages.

**Three designs**
The local method can calculate results for block, staggered, and general
designs. The triple robustness proof is only for a block design.
:::

## Know the limits

- The triple robustness result has a proof only for a block design.
- The bootstrap interval is reliable only with many treated units. fasttrop
  gives a warning message below five treated units.
- The global method works only for a block design.
- The package cannot use covariates.
- In our simulations with larger trends for the treated units, TROP reduced the
  bias, but some bias stayed. In a staggered design, the bias was 20% smaller
  than the bias of TWFE imputation.
- Our simulations with a staggered design have 60 units, four cohorts, and one
  tuning grid. We did not compare TROP with other modern estimators on real
  data.
- fasttrop follows arXiv version 3 of the paper. A line-by-line check against
  the journal version is not complete.

The article [What TROP can and cannot do](articles/limits.html) explains each
limit.

## Install the package

fasttrop is not on CRAN at this time. To install the package from GitHub, do
these steps:

1. Make sure that you have a C compiler. On Windows, Rtools has the compiler.
2. If you do not have the `remotes` package, install it.
3. Run this R code:

```r
remotes::install_github("farach/fasttrop")
```

## A first example

```r
library(fasttrop)

# A practice panel: 8 of 40 units get the treatment in the last 5 periods.
panel <- trop_simulate(n_units = 40, n_periods = 20, n_treated = 8, n_post = 5,
                       seed = 1)

fit <- trop(panel, outcome = "y", treatment = "d", unit = "unit", time = "time",
            cv_cells = 200, n_boot = 200, seed = 1)
summary(fit)
```

The article [Get started](articles/fasttrop.html) explains each step and the
output.

## Who made what

The TROP method is from Susan Athey, Guido Imbens, Zhaonan Qu, and Davide
Viviano (2026). fasttrop is an independent R package by Alex Farach. The
authors of the method did not write this package. If you use fasttrop, cite
the paper. The [citation page](authors.html#citation) gives the references.
