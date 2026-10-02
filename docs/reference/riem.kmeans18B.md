# K-Means Clustering with Lightweight Coreset

Apply weighted Lloyd iterations to an independently sampled lightweight
coreset of manifold-valued observations. The Euclidean coreset
approximation theorem is not asserted for arbitrary manifold geometries.

## Usage

``` r
riem.kmeans18B(
  riemobj,
  k = 2,
  M = max(k, ceiling(length(riemobj$data)/2)),
  geometry = c("intrinsic", "extrinsic"),
  ...
)
```

## Arguments

- riemobj:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- k:

  the number of clusters.

- M:

  integer number of independent draws, at least \\k\\ (default: the
  larger of \\k\\ and \\\lceil N/2 \rceil\\). Values greater than \\N\\
  are permitted.

- geometry:

  (case-insensitive) name or saved specification of a geometry
  supporting means. Legacy aliases `"intrinsic"` and `"extrinsic"` are
  accepted; see
  [`riem.geometry`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md).

- ...:

  extra parameters including

  maxiter

  :   maximum number of weighted Lloyd passes per start (default:50).

  nstart

  :   the number of random starts (default: 5).

## Value

a named list containing

- cluster:

  a length-\\N\\ vector of class labels (from \\1:k\\).

- means:

  a 3d array where each slice along 3rd dimension is a matrix
  representation of class mean.

- score:

  unweighted full-data within-cluster sum of squares (WCSS), used to
  select the best start.

- coreset:

  the selected start's sampled indices and importance weights.

- iterations,converged,termination:

  iteration count, convergence flag, and stopping reason for the
  selected start.

- objective_history:

  weighted coreset costs, including the initialization.

- empty_clusters:

  cluster indices unoccupied by the full data.

- starts:

  per-start validity, full-data score, convergence, and errors.

## Details

Sampling is with replacement, using the same probabilities and weights
as `riem.coreset18B`. Repeated indices remain separate weighted
observations. Each start uses weighted squared-distance initialization
and weighted cluster means. Mean solves use at most 200 iterations and
tolerance `1e-8`; failed mean solves invalidate that start. A Lloyd step
that materially increases weighted coreset cost is rejected. No global
minimum is guaranteed. Final labels assign every original observation to
its nearest stored center, with distance ties choosing the lowest
cluster index. A coreset can contain fewer than \\k\\ distinct
locations; empty centers are retained rather than redrawing and
conditioning the sampling distribution. A warning identifies an empty
cluster or a nonconverged selected start. Use
[`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## References

Bachem O, Lucic M, Krause A (2018). “Scalable k -Means Clustering via
Lightweight Coresets.” In *Proceedings of the 24th ACM SIGKDD
International Conference on Knowledge Discovery & Data Mining*,
1119–1127. ISBN 978-1-4503-5552-0.

## See also

[`riem.coreset18B`](https://www.kisungyou.com/Riemann/reference/riem.coreset18B.md)

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

## TRY DIFFERENT SIZES OF CORESET WITH K=3 FIXED
core1 = riem.kmeans18B(myriem, k=3, M=5)
core2 = riem.kmeans18B(myriem, k=3, M=10)
core3 = riem.kmeans18B(myriem, k=3, M=15)

## MDS FOR VISUALIZATION
mds2d = riem.mds(myriem, ndim=2)$embed

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
plot(mds2d, pch=19, main="true label", col=mylabs)
plot(mds2d, pch=19, main="kmeans18B: M=5",  col=core1$cluster)
plot(mds2d, pch=19, main="kmeans18B: M=10", col=core2$cluster)
plot(mds2d, pch=19, main="kmeans18B: M=15", col=core3$cluster)

par(opar)
```
