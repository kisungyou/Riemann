# Prepare Data on Sphere

The unit hypersphere (sphere, for short) is one of the most fundamental
curved space in studying geometry. Precisely, we denote \\(p-1)\\ sphere
in \\\mathbf{R}^p\\ by \$\$\mathcal{S}^{p-1} = \lbrace x \in
\mathbf{R}^p ~ \vert ~ x^\top x = \\x\\^2 = 1 \rbrace\$\$ where vectors
are of unit norm. In `wrap.sphere`, normalization is applied when each
data point is not on the unit sphere.

## Usage

``` r
wrap.sphere(input)
```

## Arguments

- input:

  data vectors to be wrapped as `riemdata` class. Following inputs are
  considered,

  matrix

  :   an \\(n \times p)\\ matrix of row observations of unit norm.

  list

  :   a length-\\n\\ list whose elements are length-\\p\\ vectors of
      unit norm.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times 1)\\ matrices in \\\mathcal{S}^{p-1}\\.

- size:

  dimension of the ambient space.

- name:

  name of the manifold of interests, *"sphere"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  Generate 5 observations in S^2 embedded in R^3.
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
test1 = wrap.sphere(d1)
test2 = wrap.sphere(d2)
```
