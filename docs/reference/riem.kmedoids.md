# K-Medoids Clustering

Given \\N\\ observations \\X_1, X_2, \ldots, X_N \in \mathcal{M}\\,
perform k-medoids clustering using pairwise distances.

## Usage

``` r
riem.kmedoids(riemobj, k = 2, geometry = c("intrinsic", "extrinsic"))
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- k:

  the number of clusters.

- geometry:

  (case-insensitive) name of geometry; either geodesic (`"intrinsic"`)
  or embedded (`"extrinsic"`) geometry.

## Value

a named list containing

- medoids:

  a length-\\k\\ vector of medoids' indices.

- cluster:

  a length-\\N\\ vector of class labels (from \\1:k\\).

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## See also

[`pam`](https://rdrr.io/pkg/cluster/man/pam.html)

## Examples

``` r
#-------------------------------------------------------------------
#          Example on Sphere : a dataset with three types
#
# class 1 : 10 perturbed data points near (1,0,0) on S^2 in R^3
# class 2 : 10 perturbed data points near (0,1,0) on S^2 in R^3
# class 3 : 10 perturbed data points near (0,0,1) on S^2 in R^3
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

## K-MEDOIDS WITH K=2,3,4
clust2 = riem.kmedoids(myriem, k=2)
clust3 = riem.kmedoids(myriem, k=3)
clust4 = riem.kmedoids(myriem, k=4)

## MDS FOR VISUALIZATION
mds2d = riem.mds(myriem, ndim=2)$embed

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
plot(mds2d, pch=19, main="true label", col=mylabs)
plot(mds2d, pch=19, main="K=2", col=clust2$cluster)
plot(mds2d, pch=19, main="K=3", col=clust3$cluster)
plot(mds2d, pch=19, main="K=4", col=clust4$cluster)

par(opar)
 
```
