# Prepare Data on Euclidean Space

Euclidean space \\\mathbf{R}^p\\ is the most common space for data
analysis, which can be considered as a Riemannian manifold with flat
metric. Since the space of matrices is isomorphic to Euclidean space
after vectorization, we consider the inputs as \\p\\-dimensional
vectors.

## Usage

``` r
wrap.euclidean(input)
```

## Arguments

- input:

  data vectors to be wrapped as `riemdata` class. Following inputs are
  considered,

  matrix

  :   an \\(n \times p)\\ matrix of row observations.

  list

  :   a length-\\n\\ list whose elements are length-\\p\\ vectors.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times 1)\\ matrices in \\\mathbf{R}^p\\.

- size:

  dimension of the ambient space.

- name:

  name of the manifold of interests, *"euclidean"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  Generate 5 observations in R^3 in Matrix and List.
#-------------------------------------------------------------------
## DATA GENERATION
d1 = array(0,c(5,3))
d2 = list()
for (i in 1:5){
  single  = stats::rnorm(3)
  d1[i,]  = single
  d2[[i]] = single
}

## RUN
test1 = wrap.euclidean(d1)
test2 = wrap.euclidean(d2)
```
