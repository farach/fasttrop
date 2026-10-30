# fasttrop: Fast Triply Robust Panel Estimators

The package calculates the Triply Robust Panel (TROP) estimator of
Athey, Imbens, Qu, and Viviano (2026). The estimator adds unit weights,
time weights, and a low-rank term with a nuclear-norm penalty to a model
with two-way fixed effects. The model fits run in C. The C code can use
more than one thread with 'OpenMP'.

## Details

Use [`trop()`](https://farach.github.io/fasttrop/reference/trop.md) for
a data frame in long format. Use
[`trop_matrix()`](https://farach.github.io/fasttrop/reference/trop_matrix.md)
for a matrix of outcomes and a matrix of treatment indicators.

## References

Athey, S., Imbens, G. W., Qu, Z., and Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) . Preprint:
arXiv:2508.21536.

## See also

Useful links:

- <https://farach.github.io/fasttrop/>

- <https://github.com/farach/fasttrop>

- Report bugs at <https://github.com/farach/fasttrop/issues>

## Author

**Maintainer**: Alex Farach <alexfarach@microsoft.com> \[copyright
holder\]

Authors:

- Alex Farach <alexfarach@microsoft.com> \[copyright holder\]
