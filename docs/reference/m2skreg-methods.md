# Methods for Scalar-Response Kernel Regression Fits

Methods for Scalar-Response Kernel Regression Fits

## Usage

``` r
# S3 method for class 'm2skreg'
fitted(object, ...)

# S3 method for class 'm2skreg'
residuals(object, ...)

# S3 method for class 'm2skreg'
summary(object, ...)

# S3 method for class 'm2skreg'
print(x, ...)

# S3 method for class 'summary.m2skreg'
print(x, ...)
```

## Arguments

- object, x:

  A fitted `m2skreg` object, or its summary for the summary printing
  method.

- ...:

  Additional arguments; currently unused.

## Value

`fitted` and `residuals` return numeric vectors; residuals are observed
responses minus training fitted values. `summary` returns a
`summary.m2skreg` list. Printing returns its input invisibly.
