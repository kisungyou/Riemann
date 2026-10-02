# Two-Sample Test modified from Biswas and Ghosh (2014)

Given \\M\\ observations \\X_1, X_2, \ldots, X_M \in \mathcal{M}\\ and
\\N\\ observations \\Y_1, Y_2, \ldots, Y_N \in \mathcal{M}\\, perform
the permutation test of equal distribution \$\$H_0~:~\mathcal{P}\_X =
\mathcal{P}\_Y\$\$ by the method from Biswas and Ghosh (2014). The
method, originally proposed for Euclidean-valued data, is adapted to the
general Riemannian manifold with intrinsic/extrinsic distance.

## Usage

``` r
riem.test2bg14(riemobj1, riemobj2, geometry = NULL, ...)
```

## Arguments

- riemobj1:

  a S3 `"riemdata"` class for \\M\\ manifold-valued data.

- riemobj2:

  a S3 `"riemdata"` class for \\N\\ manifold-valued data.

- geometry:

  A geometry name or saved specification.

- ...:

  extra parameters including

  nperm

  :   the number of permutations (default: 999).

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

The statistic compares within-group and cross-group mean distances.
Random-label calibration requires exchangeability under the
equal-distribution null, with independent observations; it is not valid
for unaccounted paired or repeated observations. Generalizing the
distance to a manifold does not establish consistency against every
alternative. Monte Carlo p-values include the observed assignment and
all ties, and cannot be zero. Each group must contain at least two
observations.

## References

Biswas M, Ghosh AK (2014). “A nonparametric two-sample test applicable
to high dimensional data.” *Journal of Multivariate Analysis*, **123**,
160–171. ISSN 0047259X.

You K, Park H (2020). “Re-visiting Riemannian geometry of symmetric
positive definite matrices for the analysis of functional connectivity.”
*NeuroImage*, 117464. ISSN 10538119.

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
#  it is expected to return a very small number.
# \donttest{
riem.test2bg14(myriem1, myriem2, nperm=999)
#> 
#>  Biswas-Ghosh distance-summary random-label test
#> 
#> data:  myriem1 and myriem2
#> Tmn = 3.6793, p-value = 0.001
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
  pvals[i] = riem.test2bg14(Xriem, Yriem, nperm=999)$p.value
}

emperr = round(sum((pvals <= 0.05))/ntest, 5)
print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
} # }
```
