# Prepare Data on Multinomial Manifold

Multinomial manifold is referred to the strictly positive data that sums
to 1. Also known as probability simplex or positive orthant, we denote
\\(p-1)\\ simplex in \\\mathbf{R}^p\\ by \$\$\Delta^{p-1} = \lbrace x
\in \mathbf{R}^p~\vert~ \sum\_{i=1}^p x_i = 1, x_i \> 0 \rbrace\$\$ in
that data are positive \\L_1\\ unit-norm vectors. In `wrap.multinomial`,
normalization is applied when each data point is not on the simplex,
accepting positive counts or masses. Zero, negative, and nonfinite
inputs are rejected before normalization; no pseudocount is added.

## Usage

``` r
wrap.multinomial(input)
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

  a list of \\(p\times 1)\\ matrices in \\\Delta^{p-1}\\.

- size:

  dimension of the ambient space.

- name:

  name of the manifold of interests, *"multinomial"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#-------------------------------------------------------------------
## DATA GENERATION
d1 = array(0,c(5,3))
d2 = list()
for (i in 1:5){
  single  = abs(stats::rnorm(3))
  d1[i,]  = single
  d2[[i]] = single
}

## RUN
test1 = wrap.multinomial(d1)
test2 = wrap.multinomial(d2)
```
