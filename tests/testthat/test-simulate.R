test_that("trop_simulate returns the documented panel and restores the RNG", {
  set.seed(5)
  before <- .Random.seed
  d <- trop_simulate(n_units = 12, n_periods = 7, n_treated = 3, n_post = 2, seed = 1)
  expect_identical(.Random.seed, before)
  expect_named(d, c("unit", "time", "y", "d", "y0"))
  expect_equal(nrow(d), 84L)
  expect_equal(sum(d$d), 6L)
  expect_equal(d$y - d$y0, d$d * 1)
  expect_identical(d, trop_simulate(n_units = 12, n_periods = 7, n_treated = 3, n_post = 2,
                                    seed = 1))
  st <- trop_simulate(n_units = 12, n_periods = 7, n_treated = 3, n_post = 3,
                      design = "staggered", seed = 2)
  expect_true(all(tapply(st$d, st$unit, function(x) all(diff(x) >= 0))))
  expect_error(trop_simulate(n_units = 5, n_treated = 5), "smaller")
})
