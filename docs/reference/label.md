# S3 method for mixture model : predict labels

Given a fitted mixture model of \\K\\ components, predict labels of
observations accordingly.

## Usage

``` r
label(object, newdata)
```

## Arguments

- object:

  a fitted mixture model of `riemmix` class.

- newdata:

  data of \\n\\ objects (vectors, matrices) that can be wrapped by one
  of `wrap.*` functions in the Riemann package.

## Value

a length-\\n\\ vector of class labels.

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
newlabel = label(k3, locations)
# }
```
