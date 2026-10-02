# Spherical Laplace Distribution

This is a collection of tools for learning with spherical Laplace (SL)
distribution on a \\(p-1)\\-dimensional sphere in \\\mathbf{R}^p\\
including sampling, density evaluation, and maximum likelihood
estimation of the parameters. The SL distribution is characterized by
the following density function, \$\$f\_{SL}(x; \mu, \sigma) =
\frac{1}{C(\sigma)} \exp \left( -\frac{d(x,\mu)}{\sigma} \right)\$\$ for
location and scale parameters \\\mu\\ and \\\sigma\\ respectively.

## Usage

``` r
dsplaplace(data, mu, sigma, log = FALSE)

rsplaplace(n, mu, sigma)

mle.splaplace(data, method = c("DE", "Optimize", "Newton"), ...)
```

## Arguments

- data:

  data vectors in form of either an \\(n\times p)\\ matrix or a
  length-\\n\\ list. See
  [`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)
  for descriptions on supported input types.

- mu:

  a length-\\p\\ unit-norm vector of location.

- sigma:

  a positive scale parameter; `Inf` gives the uniform limit.

- log:

  a logical; `TRUE` to return log-density, `FALSE` for densities without
  logarithm applied.

- n:

  the number of samples to be generated.

- method:

  an algorithm name for scale parameter estimation. It should be one of
  `"Newton"`, `"Optimize"`, and `"DE"` (case-sensitive).

- ...:

  extra parameters for computations, including

  maxiter

  :   iteration budget for each location, likelihood-search, or
      polishing stage; rounded and raised to at least 10 (default: 50).

  eps

  :   positive tolerance, capped at 1e-6; the relative likelihood-score
      tolerance is additionally floored at 1e-10 (default: 1e-6).

  use.exact

  :   for Newton, use moment-based derivatives (`TRUE`) or finite
      differences of the score (`FALSE`, default). Both use adaptive
      radial quadrature.

## Value

`dsplaplace` gives a vector of evaluated densities given samples.
`rsplaplace` generates unit-norm vectors in \\\mathbf{R}^p\\ wrapped in
a list. `mle.splaplace` computes MLEs and returns a list containing
estimates of location (`mu`) and scale (`sigma`) parameters.

## Details

Scale estimation uses an adaptively bracketed likelihood in the
reciprocal scale. Newton updates are safeguarded by that bracket;
Optimize and DE search the log reciprocal scale and use safeguarded
Newton polishing if their score has not reached the requested tolerance.
The uniform boundary is returned as `sigma = Inf`. For coincident
observations the likelihood is unbounded and `sigma = 0` is returned
with a warning; this point-mass limit has no surface-area density and is
not accepted by `dsplaplace` or `rsplaplace`. The location is obtained
by local intrinsic-median optimization and need not be globally optimal
for data spread across the sphere. Failure of the location iteration to
converge produces a warning; the returned scale then optimizes the
likelihood conditional on that last location estimate. Log densities use
the log kernel and logarithmic normalizer directly.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## Examples

``` r
# \donttest{
# -------------------------------------------------------------------
#          Example with Spherical Laplace Distribution
#
# Given a fixed set of parameters, generate samples and acquire MLEs.
# Especially, we will see the evolution of estimation accuracy.
# -------------------------------------------------------------------
## DEFAULT PARAMETERS
true.mu  = c(1,0,0,0,0)
true.sig = 1

## GENERATE A RANDOM SAMPLE OF SIZE N=1000
big.data = rsplaplace(1000, true.mu, true.sig)

## ITERATE FROM 50 TO 1000 by 10
idseq = seq(from=50, to=1000, by=10)
nseq  = length(idseq)

hist.mu  = rep(0, nseq)
hist.sig = rep(0, nseq)

for (i in 1:nseq){
  small.data = big.data[1:idseq[i]]             # data subsetting
  small.MLE  = mle.splaplace(small.data)        # compute MLE
  
  hist.mu[i]  = acos(sum(small.MLE$mu*true.mu)) # difference in mu
  hist.sig[i] = small.MLE$sigma
}

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(idseq, hist.mu,  "b", pch=19, cex=0.5, 
     main="difference in location", xlab="sample size")
plot(idseq, hist.sig, "b", pch=19, cex=0.5, 
     main="scale parameter", xlab="sample size")
abline(h=true.sig, lwd=2, col="red")

par(opar)
# }
```
