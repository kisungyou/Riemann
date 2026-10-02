# Compute Pairwise Distances for Two Sets of Data

Given \\M\\ observations \\X_1, X_2, \ldots, X_M \in \mathcal{M}\\ and
\\N\\ observations \\Y_1, Y_2, \ldots, Y_N \in \mathcal{M}\\, compute
pairwise distances between two sets' elements.

## Usage

``` r
riem.pdist2(riemobj1, riemobj2, geometry = NULL)
```

## Arguments

- riemobj1:

  a S3 `"riemdata"` class for \\M\\ manifold-valued data.

- riemobj2:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- geometry:

  A name or saved specification accepted by
  [`riem.geometry`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md).

## Value

an \\(M\times N)\\ matrix of distances.

## Examples

``` r
# \donttest{
#-------------------------------------------------------------------
#          Example on Sphere : a dataset with two types
#
#  group1 : 10 perturbed data points near (0,0,1) on S^2 in R^3
#  group2 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#           10 perturbed data points near (0,1,0) on S^2 in R^3
#           10 perturbed data points near (0,0,1) on S^2 in R^3
#-------------------------------------------------------------------
## GENERATE DATA
mydata1 = list()
mydata2 = list()
for (i in 1:10){
  tgt = c(stats::rnorm(2, sd=0.1), 1)
  mydata1[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 1:10){
  tgt = c(1, stats::rnorm(2, sd=0.1))
  mydata2[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 11:20){
  tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
  mydata2[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 21:30){
  tgt = c(stats::rnorm(2, sd=0.1), 1)
  mydata2[[i]] = tgt/sqrt(sum(tgt^2))
}
myriem1 = wrap.sphere(mydata1)
myriem2 = wrap.sphere(mydata2)

## COMPARE TWO DISTANCES
dint = riem.pdist2(myriem1, myriem2, geometry="intrinsic")
dext = riem.pdist2(myriem1, myriem2, geometry="extrinsic")

## VISUALIZE
opar = par(no.readonly=TRUE)
par(mfrow=c(1,2))
image(dint[nrow(dint):1,], main="intrinsic", axes=FALSE)
image(dext[nrow(dext):1,], main="extrinsic", axes=FALSE)

par(opar)
# }
```
