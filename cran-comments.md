## Test environments

- Local: Windows 11 Enterprise on an ARM64 processor, build 26300, with x64
  R 4.6.1 under emulation
- GitHub Actions: Ubuntu with R-devel, R-release, R-oldrel-1, and R 4.1,
  macOS with R-release, and Windows with R-release

## R CMD check results

The local `R CMD check --as-cran --no-manual` gave the status OK with 0
errors, 0 warnings, and 0 notes. The local check did not use the remote CRAN
incoming checks. On GitHub Actions, `R CMD check --as-cran --no-manual` gave
the status OK on all six platforms.

The local TeX installation does not have `makeindex` or the Courier fonts.
Thus, the local check did not make the PDF manual. With the Latin Modern
fonts and without an index, `R CMD Rd2pdf` made the manual without LaTeX
errors.

This is a new submission. Thus, the CRAN incoming check gives one NOTE for a
new submission.

## Method reference

The method is from Athey, Imbens, Qu, and Viviano (2026)
<doi:10.1002/jae.70061>. The Description field of the DESCRIPTION file has
this reference.

## Information for the reviewers

- On Windows, the C code calls `SetThreadInformation()` for its own threads.
  The call stops execution-speed power throttling (EcoQoS) while the C code
  runs. When the C code stops, Windows controls the threads again. The package
  does not change the configuration of the operating system. With
  `trop_control(throttle = TRUE)`, the package does not make this call.
- The examples and the tests use at most 2 threads. The function
  `trop_threads()` gives at most 2 threads when the environment variable
  `_R_CHECK_LIMIT_CORES_` has a value other than "false".
- The `seed` arguments are optional. The functions do not change the state of
  the random number generator of the user.
