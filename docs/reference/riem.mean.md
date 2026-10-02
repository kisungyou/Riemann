# Fréchet Mean and Variation

Given \\N\\ observations \\X_1, X_2, \ldots, X_N \in \mathcal{M}\\,
compute Fréchet mean and variation with respect to the geometry by
minimizing \$\$\textrm{min}\_x \sum\_{n=1}^N w_n \rho^2 (x, x_n),\quad
x\in\mathcal{M}\$\$ where \\\rho (x, y)\\ is a distance for two points
\\x,y\in\mathcal{M}\\. If non-uniform weights are given, normalized
version of the mean is computed and if `weight=NULL`, it automatically
sets equal weights (\\w_i = 1/n\\) for all observations.

## Usage

``` r
riem.mean(riemobj, weight = NULL, geometry = NULL, ...)
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- weight:

  Finite nonnegative observation weights with positive sum. They are
  normalized internally. `NULL` uses equal weights. Zero-weight
  observations are validated but do not enter the calculation.

- geometry:

  Geometry name or specification. `NULL` selects the default geometry.
  Legacy `"intrinsic"` and `"extrinsic"` aliases remain available for
  supported combinations. SPD names include `"affine_invariant"` and
  `"log_euclidean"`.

- ...:

  Named controls:

  maxiter

  :   Positive maximum number of accepted iterations (default 50).

  eps

  :   Positive absolute tolerance for the compatible gradient norm
      (default `1e-5`).

  max_backtrack

  :   Positive maximum number of trial steps per iteration, at most 1024
      (default 50).

  trace

  :   Whether to retain the iteration history (default `FALSE`).

  init

  :   Optional initial matrix with the observation dimensions.
      Closed-form calculations do not use an initializer.

  Unknown controls are errors.

## Value

A `riem_summary` object retaining `mean` and `variation`, plus
`objective`, resolved `geometry`, normalized `weights`, `converged`,
`termination`, `iterations`, `gradient_norm`, `step_norm`, controls, and
an optional `trace`. `variation` equals the normalized weighted sum of
squared distances at the returned mean. An unsuccessful iteration warns
and returns the last accepted estimate. For a curved extrinsic
embedding, diagnostics and trace refer to the ambient calculation, as
recorded in `diagnostic_scope`.

## Details

Intrinsic iteration uses the average logarithm direction and an Armijo
backtracking rule for the stated squared-distance objective. Euclidean
and log-Euclidean SPD means use closed forms. On spaces with nonconvex
objectives, termination with `"stationary"` establishes first-order
stationarity to the requested tolerance only. A stationary point may be
a local minimum, a saddle, or a local maximum; the criterion does not
establish local minimality, uniqueness, or global optimality.

## Examples

``` r
#-------------------------------------------------------------------
#        Example on Sphere : points near (0,1) on S^1 in R^2
#-------------------------------------------------------------------
## GENERATE DATA
ndata = 50
mydat = array(0,c(ndata,2))
for (i in 1:ndata){
  tgt = c(stats::rnorm(1, sd=2), 1)
  mydat[i,] = tgt/sqrt(sum(tgt^2))
}
myriem = wrap.sphere(mydat)

## COMPUTE TWO MEANS
mean.int = as.vector(riem.mean(myriem, geometry="intrinsic")$mean)
mean.ext = as.vector(riem.mean(myriem, geometry="extrinsic")$mean)

## VISUALIZE
opar <- par(no.readonly=TRUE)
plot(mydat[,1], mydat[,2], pch=19, xlim=c(-1.1,1.1), ylim=c(0,1.1),
     main="BLUE-extrinsic vs RED-intrinsic")
arrows(x0=0,y0=0,x1=mean.int[1],y1=mean.int[2],col="red")
arrows(x0=0,y0=0,x1=mean.ext[1],y1=mean.ext[2],col="blue")

par(opar)
```
