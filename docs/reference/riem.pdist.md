# Compute Pairwise Distances for Data

Given \\N\\ observations \\X_1, X_2, \ldots, X_N \in \mathcal{M}\\,
compute pairwise distances.

## Usage

``` r
riem.pdist(riemobj, geometry = NULL, as.dist = FALSE)
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- geometry:

  A name or saved specification accepted by
  [`riem.geometry`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md).
  NULL uses the object's saved geometry, or intrinsic geometry for an
  ordinary wrapper.

- as.dist:

  logical; if `TRUE`, it returns `dist` object, else it returns a
  symmetric matrix.

## Value

a S3 `dist` object or \\(N\times N)\\ symmetric matrix of pairwise
distances according to `as.dist` parameter.

## Examples

``` r
#-------------------------------------------------------------------
#          Example on Sphere : a dataset with two types
#
#  group1 : perturbed data points near (0,0,1) on S^2 in R^3
#  group2 : perturbed data points near (1,0,0) on S^2 in R^3
#-------------------------------------------------------------------
## GENERATE DATA
mydata = list()
sdval  = 0.1
for (i in 1:10){
  tgt = c(stats::rnorm(2, sd=sdval), 1)
  mydata[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 11:20){
  tgt = c(1, stats::rnorm(2, sd=sdval))
  mydata[[i]] = tgt/sqrt(sum(tgt^2))
}
myriem = wrap.sphere(mydata)

## COMPARE TWO DISTANCES
dint = riem.pdist(myriem, geometry="intrinsic", as.dist=FALSE)
dext = riem.pdist(myriem, geometry="extrinsic", as.dist=FALSE)

## VISUALIZE
opar = par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
image(dint[,nrow(dint):1], main="intrinsic", axes=FALSE)
image(dext[,nrow(dext):1], main="extrinsic", axes=FALSE)

par(opar)
```
