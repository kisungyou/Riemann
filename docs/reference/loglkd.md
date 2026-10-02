# S3 method for mixture model : log-likelihood

Given a fitted mixture model \\f(x)\\ and observations \\x_1, \ldots,
x_n \in \mathcal{M}\\, compute the log-likelihood \$\$L = \log
\prod\_{i=1}^n f(x_i) = \sum\_{i=1}^n \log f(x_i)\$\$.

## Usage

``` r
loglkd(object, newdata)
```

## Arguments

- object:

  a fitted mixture model of `riemmix` class.

- newdata:

  data of \\n\\ objects (vectors, matrices) that can be wrapped by one
  of `wrap.*` functions in the Riemann package.

## Value

the log-likelihood.

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
# ---------------------------------------------------- #
#            FIT A MODEL & APPLY THE METHOD
# ---------------------------------------------------- #
# Load the 'city' data and wrap as 'riemobj'
data(cities)
locations = cities$cartesian
embed2    = array(0,c(60,2)) 
for (i in 1:60){
   embed2[i,] = sphere.xyz2geo(locations[i,])
}

# Fit a model
k3 = moSN(locations, k=3)

# Evaluate
newloglkd = round(loglkd(k3, locations), 3)
print(paste0("Log-likelihood for K=3 model fit : ", newloglkd))
#> [1] "Log-likelihood for K=3 model fit : 88.582"
# }
```
