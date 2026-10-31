# Glossary

The articles on this website use these words. Each definition gives the
meaning for TROP and for fasttrop.

- ATT: The average treatment effect on the treated. It is the mean of
  the effects of the treatment on all treated cells.

- Bias: The difference between the mean of many estimates and the true
  value. An estimator with a small bias is correct on average.

- Block design: A design in which all treated units start the treatment
  in the same period and stay treated to the end of the panel.

- Bootstrap: A method that calculates the uncertainty of an estimate. It
  takes many samples of units with replacement and calculates the
  estimate again for each sample. fasttrop uses a unit bootstrap.

- Cell: One unit in one period.

- Confidence interval: A range of values for the true effect. A 95%
  interval contains the true effect in approximately 95% of samples, if
  the assumptions of the method are correct.

- Coverage: The share of simulated samples in which the confidence
  interval contains the true value. For a 95% interval, the target is
  95%.

- Cross-validation: A method that selects the tuning parameters. It
  hides data with a known outcome, calculates the hidden data from the
  other data, and measures the error. Then it selects the tuning
  parameters with the smallest error.

- Difference in differences (DiD): A method that compares the change of
  the outcome over time in the treated units with the change in the
  untreated units. It assumes parallel trends.

- Factor model: A model in which hidden patterns over time drive the
  outcomes of all units. Each unit follows each pattern with its own
  strength, its loading.

- General design: A design in which the treatment stops for one or more
  units.

- Global method: The method of the Python package from the authors of
  TROP. It calculates one weighted fit for all treated cells. It works
  only for a block design.

- Leave-one-out criterion: The cross-validation of the local method. For
  each untreated cell in the criterion, it calculates a fit without that
  cell. Then it measures the error of the fit for that cell.

- Loading: The strength with which a unit follows one hidden pattern of
  a factor model.

- Local method: The method of Eq. 12 and Algorithm 2 of the paper. Each
  treated cell has its own unit weights, time weights, and fit. It works
  for block, staggered, and general designs. It is the default method of
  fasttrop.

- Low-rank term: The part of the model that holds the hidden patterns. A
  low rank means a small number of patterns.

- Matrix completion: A method that estimates outcomes that the data does
  not show, from patterns in the observed outcomes. It uses a low-rank
  model.

- Monte Carlo simulation: A study with many simulated data sets with a
  known truth. It measures how well an estimator works.

- Monte Carlo standard error: The uncertainty of a simulation result. It
  comes from the limited number of simulated data sets.

- Normal interval: The estimate plus and minus a normal quantile times
  the standard error. For a 95% interval, the quantile is 1.96. This is
  the default interval of fasttrop.

- Nuclear-norm penalty:

  The penalty that controls the low-rank term. A larger penalty gives
  fewer hidden patterns. With an infinite penalty (`nn = Inf`), the
  model has no low-rank term.

- Outcome without the treatment: The outcome of a treated cell if the
  unit did not get the treatment. The data cannot show this outcome.
  TROP estimates it from the untreated cells.

- Panel data: Data with an outcome for many units in many periods.

- Parallel trends: The assumption that, without the treatment, the
  treated units and the untreated units change in the same way over
  time.

- Percentile interval: An interval from the quantiles of the bootstrap
  estimates. A 95% interval uses the 2.5% and 97.5% quantiles.

- Placebo cross-validation: The cross-validation of the global method.
  It gives a false treatment to some untreated units and estimates its
  effect. The true effect of a false treatment is zero, so the estimate
  shows the error.

- RMSE: The root mean square error. It is the square root of the mean of
  the squared differences between the estimates and the true value. It
  includes both the bias and the noise of the estimates.

- Staggered design: A design in which the treated units start the
  treatment in different periods and stay treated.

- Standard error: A measure of the uncertainty of an estimate. A smaller
  standard error shows a more precise estimate.

- Synthetic control: A method that makes a weighted mix of untreated
  units. The mix copies the history of the treated unit before the
  treatment.

- Synthetic difference in differences (SDID): A method that adds unit
  weights and time weights to difference in differences.

- Thread: One path of calculation on a computer. With more threads,
  fasttrop can calculate more fits at the same time.

- Time weights:

  Weights that give more importance to periods near the treated period.
  The tuning parameter `time` controls how fast they decrease.

- Treated cell: A cell in which the unit has the treatment.

- Triple robustness: The property of TROP that its bias can be small
  when only one of its three parts works well. The three parts are the
  unit weights, the time weights, and the model. The paper proves this
  property for a block design with fixed weights, under its assumptions.

- Tuning parameters:

  The three numbers that control TROP: `time`, `unit`, and `nn`.
  Cross-validation selects them, or you can set them.

- TWFE imputation: A method with two steps. First, it fits an effect for
  each unit and an effect for each period to the untreated cells only.
  Second, it uses this model to calculate the outcome without the
  treatment for each treated cell. TROP with equal weights and no
  low-rank term gives the same result.

- Two-way fixed effects (TWFE): A model with an effect for each unit and
  an effect for each period. A conventional TWFE regression fits this
  model to all cells, with one coefficient for the treatment. In a
  staggered design, its result can be different from the result of TWFE
  imputation.

- Unit: One member of the panel, for example a state, a firm, or a
  person.

- Unit weights:

  Weights that give more importance to units with an outcome history
  that is similar to the treated unit. The tuning parameter `unit`
  controls how fast they decrease.

- Untreated cell: A cell in which the unit does not have the treatment.
  Before the treatment starts, the cells of a treated unit are also
  untreated cells.
