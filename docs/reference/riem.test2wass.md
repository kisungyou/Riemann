# Two-Sample Test with Wasserstein Metric

Given \\M\\ observations \\X_1, X_2, \ldots, X_M \in \mathcal{M}\\ and
\\N\\ observations \\Y_1, Y_2, \ldots, Y_N \in \mathcal{M}\\,
permutation test based on the Wasserstein metric (see
[`riem.wasserstein`](https://www.kisungyou.com/Riemann/reference/riem.wasserstein.md)
for more details) is applied to test whether two distributions are same
or not, i.e., \$\$H_0~:~\mathcal{P}\_X = \mathcal{P}\_Y\$\$ with
Wasserstein metric \\\mathcal{W}\_p\\ being the measure of discrepancy
between two samples.

## Usage

``` r
riem.test2wass(riemobj1, riemobj2, p = 2, geometry = NULL, ...)
```

## Arguments

- riemobj1:

  a S3 `"riemdata"` class for \\M\\ manifold-valued data.

- riemobj2:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- p:

  an exponent for Wasserstein distance \\\mathcal{W}\_p\\ (default: 2).

- geometry:

  (case-insensitive) name of geometry; either geodesic (`"intrinsic"`)
  or embedded (`"extrinsic"`) geometry.

- ...:

  extra parameters including

  nperm

  :   the number of permutations (default: 999).

  use.smooth

  :   a logical, default `FALSE`. `TRUE` requests an experimental IPOT
      approximation, without an optimality certificate.

## Value

a (list) object of `S3` class `htest` containing:

- statistic:

  a test statistic.

- p.value:

  \\p\\-value under \\H_0\\.

- alternative:

  alternative hypothesis.

- method:

  name of the test.

- data.name:

  name(s) of provided sample data.

## Details

Random label permutations require independent observations that are
exchangeable under the identical-distributions null. This procedure does
not handle paired, clustered or repeated observations. The Monte Carlo
p-value includes ties and uses the plus-one correction. All assignments
use the same solver and exponent. Saved geometry, permuted statistics,
and the Monte Carlo standard-error diagnostic are retained.

## Examples

``` r
#-------------------------------------------------------------------
#          Example on Sphere : a dataset with two types
#
# class 1 : 20 perturbed data points near (1,0,0) on S^2 in R^3
# class 2 : 30 perturbed data points near (0,1,0) on S^2 in R^3
#-------------------------------------------------------------------
## GENERATE DATA
mydata1 = list()
mydata2 = list()
for (i in 1:20){
  tgt = c(1, stats::rnorm(2, sd=0.1))
  mydata1[[i]] = tgt/sqrt(sum(tgt^2))
}
for (i in 1:20){
  tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
  mydata2[[i]] = tgt/sqrt(sum(tgt^2))
}
myriem1 = wrap.sphere(mydata1)
myriem2 = wrap.sphere(mydata2)

## PERFORM PERMUTATION TEST
#  it is expected to return a very small number, but 
#  small number of 'nperm' may not give a reasonable p-value.
# \donttest{
riem.test2wass(myriem1, myriem2, nperm=99, use.smooth=FALSE)
#> 
#>  Wasserstein two-sample permutation test
#> 
#> data:  Supplied samples
#> Wmn = 1.5732, p-value = 0.01
#> alternative hypothesis: true exchangeable_group_distributions is  0
#> 
# }

if (FALSE) { # \dontrun{
## CHECK WITH EMPIRICAL TYPE-1 ERROR
set.seed(777)
ntest = 1000
pvals = rep(0,ntest)

for (i in 1:ntest){
  X = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
  Y = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
  Xnorm = X/sqrt(rowSums(X^2))
  Ynorm = Y/sqrt(rowSums(Y^2))
  
  Xriem = wrap.sphere(Xnorm)
  Yriem = wrap.sphere(Ynorm)
  pvals[i] = riem.test2wass(Xriem, Yriem, nperm=999)$p.value
  print(paste0("iteration ",i,"/",ntest," complete.."))
}

emperr = round(sum((pvals <= 0.05))/ntest, 5)
print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
} # }
```
