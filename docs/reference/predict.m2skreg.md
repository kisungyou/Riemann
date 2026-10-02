# Prediction for Manifold-to-Scalar Kernel Regression

Predicts under the fitted geometry and bandwidth without refitting or
using prediction-batch statistics. Training predictors and responses are
required.

## Usage

``` r
# S3 method for class 'm2skreg'
predict(
  object,
  newdata,
  geometry = NULL,
  diagnostics = FALSE,
  block_size = 256L,
  ...
)
```

## Arguments

- object:

  A fitted `m2skreg` object.

- newdata:

  A compatible `riemdata` object containing new predictors.

- geometry:

  Usually `NULL`, which uses the fitted geometry. An explicit value must
  resolve to the same geometry. An old serialized object without
  geometry metadata requires an explicit value and produces a warning,
  because its original geometry cannot be inferred reliably.

- diagnostics:

  Whether to return support diagnostics with predictions.

- block_size:

  Positive integer limiting the number of new observations in each
  cross-distance calculation. It does not change the fitted model.

- ...:

  Reserved for future arguments; unknown arguments are rejected.

## Value

By default, a numeric vector with one prediction per observation. With
`diagnostics=TRUE`, a list containing `prediction` and a data frame
`diagnostics` with `effective_n`, `nearest_distance`, and
`nearest_over_bandwidth`. The last quantity may be infinite when the
ratio is not representable even though the prediction remains finite.

## See also

[`riem.m2skreg`](https://www.kisungyou.com/Riemann/reference/riem.m2skreg.md)
