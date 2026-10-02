# Frechet Analysis of Variance

Compares population Frechet means and variances using the statistic of
Dubey and Muller (2019), equation (14). The asymptotic chi-squared
calibration targets equality of those summaries, not arbitrary equality
of distributions. The permutation version requires exchangeable
observations under the stronger null of identical group distributions,
and targets alternatives detected by the same mean-and-variance
statistic.

## Usage

``` r
riem.fanova(..., maxiter = 50, eps = 1e-05, geometry = NULL)

riem.fanovaP(..., maxiter = 50, eps = 1e-05, nperm = 99, geometry = NULL)
```

## Arguments

- ...:

  At least two compatible `riemdata` objects, each containing at least
  three observations. Observations must be independent within and across
  groups for the documented calibration; paired, repeated or clustered
  observations require a separate restricted resampling procedure.

- maxiter:

  Positive integer iteration budget for each Frechet mean.

- eps:

  Finite positive tolerance for the shared mean solver.

- geometry:

  A geometry name or saved specification with a compatible mean and
  distance implementation. By default all groups must resolve to the
  same geometry.

- nperm:

  Positive integer number of random label permutations.

## Value

An `htest` object with statistic, p-value, null/alternative, resolved
geometry, calibration, group sizes and mean diagnostics. The permutation
version additionally retains permuted statistics, `nperm`, and `mc_se`.

## Details

The asymptotic result requires unique population/sample means, positive
variances of the squared distances to each group mean, suitable
moment/metric-entropy conditions, and group proportions bounded away
from zero. The cited paper gives sufficient bounded-space conditions;
accepting an input manifold does not establish them for every population
distribution. The function rejects failed mean fits and a degenerate
empirical denominator.

For group proportions \\\lambda_j=n_j/n\\, the two denominators are
\\\sum_j\lambda_j/\widehat\sigma_j^2\\ and
\\\sum_j\lambda_j^2\widehat\sigma_j^2\\. Distances are jointly rescaled
before forming the statistic; this leaves it unchanged and avoids
unit-scale overflow. Numerical rescaling does not transform the
underlying geometry.

`riem.fanovaP` refits all group means for every sampled label
assignment. It uses \\(1+\\\\T_b\geq T\_{obs}\\)/(nperm+1)\\, including
ties. The pooled fit is unchanged by label permutation and is reused.
Permutations follow R's current random-number state. The returned Monte
Carlo standard error is a plug-in diagnostic for simulation variability,
not inferential uncertainty in the scientific effect.

## References

Dubey P, Müller H (2019). “Fréchet analysis of variance for random
objects.” *Biometrika*, **106**(4), 803–821. ISSN 0006-3444, 1464-3510.

## Examples

``` r
X <- wrap.euclidean(matrix(c(-1, 0, 2, 3), ncol = 1))
Y <- wrap.euclidean(matrix(c(0, 1, 4, 8), ncol = 1))
riem.fanova(X, Y)
#> 
#>  Frechet Analysis of Variance on euclidean
#> 
#> data:  Supplied groups
#> Tn = 3.6818, df = 1, p-value = 0.05501
#> alternative hypothesis: true equal_frechet_means_and_variances is  0
#> 
set.seed(17)
riem.fanovaP(X, Y, nperm = 19)
#> 
#>  Frechet Analysis of Variance on euclidean (permutation calibration)
#> 
#> data:  Supplied groups
#> Tn = 3.6818, df = 1, p-value = 0.35
#> alternative hypothesis: true exchangeable_group_distributions is  0
#> 
```
