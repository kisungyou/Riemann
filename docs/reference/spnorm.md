# Spherical Normal Distribution

We provide tools for an isotropic spherical normal (SN) distributions on
a \\(p-1)\\-sphere in \\\mathbf{R}^p\\ for sampling, density evaluation,
and maximum likelihood estimation of the parameters where the density is
defined as \$\$f\_{SN}(x; \mu, \lambda) = \frac{1}{Z(\lambda)} \exp
\left( -\frac{\lambda}{2} d^2(x,\mu) \right)\$\$ for location and
concentration parameters \\\mu\\ and \\\lambda\\ respectively and the
normalizing constant \\Z(\lambda)\\.

## Usage

``` r
dspnorm(data, mu, lambda, log = FALSE)

rspnorm(n, mu, lambda)

mle.spnorm(data, method = c("Newton", "Halley", "Optimize", "DE"), ...)
```

## Arguments

- data:

  data vectors in form of either an \\(n\times p)\\ matrix or a
  length-\\n\\ list. See
  [`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)
  for descriptions on supported input types.

- mu:

  a length-\\p\\ unit-norm vector of location.

- lambda:

  a finite nonnegative concentration parameter; zero is the uniform
  distribution.

- log:

  a logical; `TRUE` to return log-density, `FALSE` for densities without
  logarithm applied.

- n:

  the number of samples to be generated.

- method:

  an algorithm name for concentration parameter estimation. It should be
  one of `"Newton"`,`"Halley"`,`"Optimize"`, and `"DE"` (case
  sensitive).

- ...:

  extra parameters for computations, including

  maxiter

  :   iteration budget for each location, likelihood-search, or
      polishing stage; rounded and raised to at least 10 (default: 50).

  eps

  :   positive tolerance, capped at 1e-6; the relative likelihood-score
      tolerance is additionally floored at 1e-10 (default: 1e-6).

## Value

`dspnorm` gives a vector of evaluated densities given samples. `rspnorm`
generates unit-norm vectors in \\\mathbf{R}^p\\ wrapped in a list.
`mle.spnorm` computes MLEs and returns a list containing estimates of
location (`mu`) and concentration (`lambda`) parameters.

## Details

Concentration estimation uses an adaptively bracketed likelihood score.
Newton and Halley updates are safeguarded by that bracket; Optimize and
DE search the log concentration and use safeguarded Newton polishing if
their score has not reached the requested tolerance. All methods check
the uniform boundary, returning `lambda = 0` when appropriate. For
coincident observations the likelihood is unbounded and `lambda = Inf`
is returned with a warning; this point-mass limit has no density with
respect to spherical surface area and is not accepted by `dspnorm` or
`rspnorm`. The location is obtained by local intrinsic-mean
optimization; a globally optimal location is not guaranteed for data
spread across the sphere. Failure of the location iteration to converge
produces a warning; the returned concentration then optimizes the
likelihood conditional on that last location estimate. Log densities are
evaluated directly, including a logarithmic normalizer.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## References

Hauberg S (2018). “Directional Statistics with the Spherical Normal
Distribution.” In *2018 21st International Conference on Information
Fusion (FUSION)*, 704–711. ISBN 978-0-9964527-6-2.

You K, Suh C (2022). “Parameter Estimation and Model-Based Clustering
with Spherical Normal Distribution on the Unit Hypersphere.”
*Computational Statistics and Data Analysis*, 107457. ISSN 01679473.

## Examples

``` r
# \donttest{
# -------------------------------------------------------------------
#          Example with Spherical Normal Distribution
#
# Given a fixed set of parameters, generate samples and acquire MLEs.
# Especially, we will see the evolution of estimation accuracy.
# -------------------------------------------------------------------
## DEFAULT PARAMETERS
true.mu  = c(1,0,0,0,0)
true.lbd = 5

## GENERATE DATA N=1000
big.data = rspnorm(1000, true.mu, true.lbd)

## ITERATE FROM 50 TO 1000 by 10
idseq = seq(from=50, to=1000, by=10)
nseq  = length(idseq)

hist.mu  = rep(0, nseq)
hist.lbd = rep(0, nseq)

for (i in 1:nseq){
  small.data = big.data[1:idseq[i]]          # data subsetting
  small.MLE  = mle.spnorm(small.data) # compute MLE
  
  hist.mu[i]  = acos(sum(small.MLE$mu*true.mu)) # difference in mu
  hist.lbd[i] = small.MLE$lambda
}

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(idseq, hist.mu,  "b", pch=19, cex=0.5, main="difference in location")
plot(idseq, hist.lbd, "b", pch=19, cex=0.5, main="concentration param")
abline(h=true.lbd, lwd=2, col="red")

par(opar)
# }
```
