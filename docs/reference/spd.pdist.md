# Pairwise Distance on SPD Manifold

Given \\N\\ observations \\X_1, X_2, \ldots, X_N\\ in SPD manifold,
compute pairwise distances among observations. Stein distances use a
scaled relative spectrum, and Wasserstein distances use an aligned
Cholesky-factor residual, avoiding determinant overflow and cancellation
of nearly equal traces. Numerically singular inputs can still prevent
these factorizations and are reported as errors.

## Usage

``` r
spd.pdist(spdobj, geometry, as.dist = FALSE)
```

## Arguments

- spdobj:

  a S3 `"riemdata"` class of SPD-valued data.

- geometry:

  name of the geometry to be used. See
  [`spd.geometry`](https://www.kisungyou.com/Riemann/reference/spd.geometry.md)
  for supported geometries.

- as.dist:

  logical; if `TRUE`, it returns a `dist` object. Else, it returns a
  symmetric matrix.

## Value

a S3 `dist` object or \\(N\times N)\\ symmetric matrix of pairwise
distances according to `as.dist` parameter.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## Examples

``` r
# \donttest{
#-------------------------------------------------------------------
#                   Two Types of Covariances
#
#  group1 : perturbed from data by N(0,1) in R^3
#  group2 : perturbed from data by [sin(x); cos(x); sin(x)*cos(x)]
#-------------------------------------------------------------------
## GENERATE DATA
spd_mats = array(0,c(3,3,20))
for (i in 1:10){
  spd_mats[,,i] = stats::cov(matrix(rnorm(50*3), ncol=3))
}
for (j in 11:20){
  randvec = stats::rnorm(50, sd=3)
  randmat = cbind(sin(randvec), cos(randvec), sin(randvec)*cos(randvec))
  spd_mats[,,j] = stats::cov(randmat + matrix(rnorm(50*3, sd=0.1), ncol=3))
}

## WRAP IT AS SPD OBJECT
spd_obj = wrap.spd(spd_mats)

## COMPUTE PAIRWISE DISTANCES
#  Geometries are case-insensitive.
pdA = spd.pdist(spd_obj, "airM")
pdL = spd.pdist(spd_obj, "lErm")
pdJ = spd.pdist(spd_obj, "Jeffrey")
pdS = spd.pdist(spd_obj, "stEin")
pdW = spd.pdist(spd_obj, "wasserstein")

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,3), pty="s")
image(pdA, axes=FALSE, main="AIRM")
image(pdL, axes=FALSE, main="LERM")
image(pdJ, axes=FALSE, main="Jeffrey")
image(pdS, axes=FALSE, main="Stein")
image(pdW, axes=FALSE, main="Wasserstein")
par(opar)

# }
```
