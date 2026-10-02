# Tangent Principal Component Analysis

Given \\N\\ observations \\X_1, X_2, \ldots, X_N \in \mathcal{M}\\,
Principal Geodesic Analysis (PGA) finds a low-dimensional embedding by
decomposing 2nd-order information in tangent space at an intrinsic mean
of the data.

## Usage

``` r
riem.pga(
  riemobj,
  ndim = 2,
  geometry = NULL,
  weight = NULL,
  maxiter = 200,
  eps = 1e-08,
  rank.tol = sqrt(.Machine$double.eps),
  center.tangent = TRUE
)

# S3 method for class 'riem_pga'
predict(object, newdata, ...)

# S3 method for class 'riem_pga'
riem.reconstruct(object, scores = object$embed, ...)

# S3 method for class 'riem_pga'
print(x, ...)

# S3 method for class 'riem_pga'
summary(object, ...)

# S3 method for class 'summary.riem_pga'
print(x, ...)

# S3 method for class 'riem_pga'
plot(x, components = c(1L, 2L), ...)
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- ndim:

  positive integer requested number of components. Requests above
  numerical rank are truncated with a warning.

- geometry:

  a supported geometry name or saved geometry specification. Defaults to
  the geometry of the data. Supports affine-invariant and log-Euclidean
  SPD, round spheres, Euclidean observations, and regular
  orthogonal-quotient landmark shapes.

- weight:

  nonnegative observation weights, normalized to sum to one.

- maxiter:

  maximum iterations for the reference mean.

- eps:

  convergence tolerance for the reference mean.

- rank.tol:

  relative singular-value tolerance for numerical rank.

- center.tangent:

  whether to subtract the weighted tangent-coordinate mean before
  decomposition. The fitted offset is retained for reuse.

- object:

  a fitted tangent PCA object.

- newdata:

  a compatible `riemdata` object.

- ...:

  additional arguments; currently unused.

- scores:

  numeric matrix with one column per retained component; by default, the
  fitted training scores are reconstructed.

- x:

  a fitted tangent PCA object.

- components:

  two retained component indices to display.

## Value

a named list containing

- center:

  an intrinsic mean in a matrix representation form.

- embed:

  an \\N\\-by-effective-rank matrix of component scores.

- loadings:

  orthonormal directions in the stored tangent coordinates.

- geometry:

  the persistent geometry specification.

- offset:

  the fitted tangent-coordinate offset.

- rank:

  numerical rank of the weighted coordinate matrix.

- variance:

  component variances in decreasing order.

- diagnostics:

  reference-mean convergence and approximation information.

## Details

This is tangent PCA, a local linear approximation, rather than an
optimization over principal geodesic submanifolds. SPD coordinates
preserve the selected metric using symmetric vectorization with
square-root-of-two off-diagonal weights. Affine-invariant coordinates
are whitened at the mean; log-Euclidean coordinates use the matrix-log
chart. Sphere coordinates use an orthonormal tangent basis. Landmark
coordinates use horizontal log vectors after alignment to the trained
reference under the full orthogonal group (reflections included);
singular alignments are rejected.

Component variances use normalized weights divided by
`1 - sum(weight^2)`, agreeing with sample PCA for equal weights. This
denominator is evaluated as twice the sum of positive pair products to
avoid cancellation with highly unequal weights. Positive weights that
underflow to zero during normalization are rejected. With only one
positive weight, the centered fit has zero rank and zero variance by
convention; this is not an estimate of a sample covariance. An
uncentered fit requires at least two positive weights. Explained
variance describes the tangent representation, not exact global manifold
variation. A nonconverged reference mean causes an error. Prediction
retains the trained mean, coordinate convention, offset and loadings.
Reconstruction returns wrapped observations in this same frame. Sphere
reconstruction is restricted to tangent norms below pi; landmark
reconstruction is restricted below pi/2 and to regular shapes.

## References

Fletcher PT, Lu C, Pizer SM, Joshi S (2004). “Principal Geodesic
Analysis for the Study of Nonlinear Statistics of Shape.” *IEEE
Transactions on Medical Imaging*, **23**(8), 995–1005. ISSN 0278-0062.

## Examples

``` r
#-------------------------------------------------------------------
#          Example on Sphere : a dataset with three types
#
# 10 perturbed data points near (1,0,0) on S^2 in R^3
# 10 perturbed data points near (0,1,0) on S^2 in R^3
# 10 perturbed data points near (0,0,1) on S^2 in R^3
#-------------------------------------------------------------------
## GENERATE DATA
mydata = list()
for (i in 1:10){
  tgt = c(1, stats::rnorm(2, sd=0.1))
  mydata[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 11:20){
  tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
  mydata[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 21:30){
  tgt = c(stats::rnorm(2, sd=0.1), 1)
  mydata[[i]] = tgt/sqrt(sum(tgt^2))
}
myriem = wrap.sphere(mydata)
mylabs = rep(c(1,2,3), each=10)

## EMBEDDING WITH MDS AND PGA
embed2mds = riem.mds(myriem, ndim=2, geometry="intrinsic")$embed
embed2pga = riem.pga(myriem, ndim=2)$embed

## VISUALIZE
opar = par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
plot(embed2mds, main="Multidimensional Scaling",    col=mylabs, pch=19)
plot(embed2pga, main="Principal Geodesic Analysis", col=mylabs, pch=19)

par(opar)
```
