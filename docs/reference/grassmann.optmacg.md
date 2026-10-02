# Estimation of Distribution Algorithm with MACG Distribution

For a function \\f : Gr(k,p) \rightarrow \mathbf{R}\\, find the
minimizer and the attained minimum value with estimation of distribution
algorithm using MACG distribution.

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

  extra parameters including

  n.start

  :   number of runs; algorithm is executed `n.start` times (default:
      10).

  maxiter

  :   maximum number of iterations for each run (default: 100).

  popsize

  :   the number of samples generated at each step for stochastic search
      (default: 100).

  ratio

  :   ratio in \\(0,1)\\ where top `ratio*popsize` samples are chosen
      for parameter update (default: 0.25).

  print.progress

  :   a logical; if `TRUE`, it prints each iteration (default: `FALSE`).

## Value

a named list containing:

- cost:

  smallest function value attained during the iterative search runs; a
  global minimum is not guaranteed.

- solution:

  a \\(p\times k)\\ matrix that attains the `cost`.

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
A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
myfunc <- function(p){                # cost function to minimize
  return(sum(-diag(t(p)%*%A%*%p)))
} 

## SOLVE THE OPTIMIZATION PROBLEM
Aout = grassmann.optmacg(myfunc, p=5, k=3, popsize=100, n.start=30)
#> Warning: The angular shape iteration did not converge; the last positive-definite estimate is returned.
#> Error: Shape update must be a finite symmetric positive-definite matrix.

## COMPUTE EIGENVALUES
#  1. USE SOLUTIONS TO THE ABOVE OPTIMIZATION 
abase   = Aout$solution
#> Error: object 'Aout' not found
eig3sol = sort(diag(t(abase)%*%A%*%abase), decreasing=TRUE)
#> Error: object 'abase' not found

#  2. USE BASIC 'EIGEN' FUNCTION
eig3dec = sort(eigen(A)$values, decreasing=TRUE)[1:3]

## VISUALIZE
opar <- par(no.readonly=TRUE)
yran = c(min(min(eig3sol),min(eig3dec))*0.95,
         max(max(eig3sol),max(eig3dec))*1.05)
#> Error: object 'eig3sol' not found
plot(1:3, eig3sol, type="b", col="red",  pch=19, ylim=yran,
     xlab="index", ylab="eigenvalue", main="compare top 3 eigenvalues")
#> Error: object 'eig3sol' not found
lines(1:3, eig3dec, type="b", col="blue", pch=19)
#> Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
legend(1.55, max(yran), legend=c("optimization","decomposition"), col=c("red","blue"),
       lty=rep(1,2), pch=19)
#> Error: object 'yran' not found
par(opar)
# }
```
