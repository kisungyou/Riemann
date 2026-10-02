# Manifold-to-Scalar Kernel Regression

Fits the Nadaraya–Watson smoother for manifold-valued predictors and
finite real scalar responses. For bandwidth \\h\>0\\, weights at \\x\\
are proportional to \\\exp\\-d(x,X_i)^2/(2h^2)\\\\. The selected
distance and training observations are retained for prediction. Training
fitted values include the observation's own weight; they are not
cross-validated predictions.

## Usage

``` r
riem.m2skreg(riemobj, y, bandwidth = 0.5, geometry = NULL)
```

## Arguments

- riemobj:

  A `riemdata` object containing the training predictors.

- y:

  A finite numeric response vector with one entry per predictor.

- bandwidth:

  One finite, strictly positive bandwidth.

- geometry:

  Geometry name or specification. `NULL` resolves the input's geometry,
  defaulting to its intrinsic geometry when unspecified. Legacy
  `"intrinsic"` and `"extrinsic"` aliases are accepted.

## Value

An object of class `m2skreg`. Legacy fields `ypred`, `bandwidth`, and
`inputs` are retained. Additional fields include resolved `geometry`,
`call`, input dimensions, schema and package versions, and
`training_diagnostics`. No distance matrix is retained.

## Details

Weights are normalized after subtracting the smallest squared distance
in the exponent. The difference is evaluated without first squaring the
distances, avoiding all-zero weights for small bandwidths or distant
predictions. Numerically negligible relative weights can still underflow
to zero. This is floating-point evaluation of the Gaussian smoother, not
a change to a nearest-neighbor model.

The effective weight count and nearest training distance divided by the
bandwidth describe the weights and distance scale. They are not
confidence intervals or assurances of adequate statistical support.

## Examples

``` r
theta <- seq(0, pi / 2, length.out = 8)
X <- wrap.sphere(cbind(cos(theta), sin(theta)))
fit <- riem.m2skreg(X, sin(2 * theta), bandwidth = 0.3)
predict(fit, X)
#> [1] 0.3097147 0.4882286 0.6854487 0.8194095 0.8194095 0.6854487 0.4882286
#> [8] 0.3097147
summary(fit)
#> Manifold-to-scalar Gaussian kernel regression
#> Observations: 8   Bandwidth: 0.3 
#> Geometry: round 
#> Training RMSE (includes self-weights): 0.1819012 
```
