# Estimation of Distribution Algorithm with MACG Distribution

Search for a subspace minimizing a basis-invariant objective by adapting
a matrix angular central Gaussian (MACG) sampling distribution to elite
frames.

## Usage

``` r
grassmann.optmacg(func, p, k, ...)
```

## Arguments

- func:

  a function to be *minimized*, returning one finite numeric value.

- p:

  dimension parameter as in \\Gr(k,p)\\.

- k:

  dimension parameter as in \\Gr(k,p)\\.

- ...:

  named controls including

  n.start

  : number of runs; algorithm is executed `n.start` times (default: 10).

  maxiter

  : maximum number of iterations for each run (default: 100).

  popsize

  : the number of samples generated at each step for stochastic search
    (default: 100).

  ratio

  : ratio in \\(0,1)\\ where top `ratio*popsize` samples are chosen for
    parameter update (default: 0.25).

  print.progress

  : a logical; if `TRUE`, it prints each iteration (default: `FALSE`).

## Value

a named list containing:

- cost:

  smallest function value attained during the iterative search runs; a
  global minimum is not guaranteed.

- solution:

  a \\(p\times k)\\ matrix that attains the `cost`.

## Details

The first run starts with uniform samples. Each additional run begins
with a uniform pilot population twice the requested size. At each
update, the best frames determine a trace-normalized MACG proposal via
at most 50 fixed-point steps. Each step mixes in a fraction
`sqrt(.Machine$double.eps)` of the identity to keep the search
distribution positive definite even when the elite subspaces coincide.
This regularization applies only to the optimizer's sampling proposal. A
run ends after five generations without improvement or at `maxiter`. For
\\k=p\\, the Grassmannian contains only one subspace.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## Examples

``` r
#-------------------------------------------------------------------
#               Optimization for Eigen-Decomposition
#
# Given (5x5) covariance matrix S, eigendecomposition is can be 
# considered as an optimization on Grassmann manifold. Here, 
# we are trying to find top 3 eigenvalues and compare.
#-------------------------------------------------------------------
# \donttest{
## PREPARE
set.seed(17)
A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
myfunc <- function(p){                # cost function to minimize
  return(sum(-diag(t(p)%*%A%*%p)))
} 

## SOLVE THE OPTIMIZATION PROBLEM
Aout = grassmann.optmacg(myfunc, p=5, k=3, popsize=100, n.start=3, maxiter=30)

## COMPUTE EIGENVALUES
#  1. USE SOLUTIONS TO THE ABOVE OPTIMIZATION 
abase   = Aout$solution
eig3sol = eigen(crossprod(abase, A %*% abase), symmetric=TRUE)$values

#  2. USE BASIC 'EIGEN' FUNCTION
eig3dec = sort(eigen(A)$values, decreasing=TRUE)[1:3]

## VISUALIZE
opar <- par(no.readonly=TRUE)
yran = c(min(min(eig3sol),min(eig3dec))*0.95,
         max(max(eig3sol),max(eig3dec))*1.05)
plot(1:3, eig3sol, type="b", col="red",  pch=19, ylim=yran,
     xlab="index", ylab="eigenvalue", main="compare top 3 eigenvalues")
lines(1:3, eig3dec, type="b", col="blue", pch=19)
legend(1.55, max(yran), legend=c("optimization","decomposition"), col=c("red","blue"),
       lty=rep(1,2), pch=19)

par(opar)
# }
```
