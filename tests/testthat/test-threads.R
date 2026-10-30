test_that("trop_threads() uses at most 2 threads under R CMD check", {
  old_opt <- options(fasttrop.threads = NULL)
  old_env <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  on.exit({
    options(old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_env)
    }
  }, add = TRUE)

  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_lte(trop_threads(), 2L)
  expect_gte(trop_threads(), 1L)

  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  expect_gte(trop_threads(), 1L)

  # The option has priority over the limit.
  options(fasttrop.threads = 3L)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_identical(trop_threads(), 3L)
})

test_that("the default control object uses trop_threads()", {
  old_opt <- options(fasttrop.threads = 1L)
  on.exit(options(old_opt), add = TRUE)
  expect_identical(trop_control()$threads, 1L)
})
