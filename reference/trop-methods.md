# Methods for `trop` objects

These methods show the results of a
[`trop()`](https://farach.github.io/fasttrop/reference/trop.md) fit and
get parts of the results.

## Usage

``` r
# S3 method for class 'trop'
print(x, digits = 4L, ...)

# S3 method for class 'trop'
summary(object, ...)

# S3 method for class 'summary.trop'
print(x, digits = 4L, ...)

# S3 method for class 'trop'
coef(object, ...)

# S3 method for class 'trop'
confint(object, parm, level = object$conf.level, ...)

# S3 method for class 'trop'
vcov(object, ...)

# S3 method for class 'trop'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x, object:

  A `trop` object.

- digits:

  The number of significant digits to show.

- ...:

  Not used.

- parm:

  Not used. The model has only one parameter, the ATT.

- level:

  The confidence level. Only the level of the fit is available, because
  the fit contains the interval.

- row.names, optional:

  Not used.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the ATT with the
name `ATT`. [`confint()`](https://rdrr.io/r/stats/confint.html) returns
a matrix with 1 row and 2 columns.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns a matrix with 1
row and 1 column.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with one row for each treated cell and its effect.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly,
and [`summary()`](https://rdrr.io/r/base/summary.html) returns a summary
object.
