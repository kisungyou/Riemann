# Prepare Data on Symmetric Positive-Definite (SPD) Manifold

The collection of symmetric positive-definite matrices is a well-known
example of matrix manifold. It is defined as \$\$\mathcal{S}\_{++}^p =
\lbrace X \in \mathbf{R}^{p\times p} ~\vert~ X^\top = X,~ v^\top
Xv\>0~\textrm{for all}~v\ne0 \rbrace\$\$ where positivity is verified by
Cholesky factorization after a scale-relative symmetry check. No ridge
or positive-definite projection is applied. Tiny physical units alone do
not make a matrix ill-conditioned. Please note that the geometry
involving semi-definite matrices is considered in `wrap.spdk`.

## Usage

``` r
wrap.spd(input)
```

## Arguments

- input:

  SPD data matrices to be wrapped as `riemdata` class. Following inputs
  are considered,

  array

  :   an \\(p\times p\times n)\\ array where each slice along 3rd
      dimension is a SPD matrix.

  list

  :   a length-\\n\\ list whose elements are \\(p\times p)\\ SPD
      matrices.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times p)\\ SPD matrices.

- size:

  size of each SPD matrix.

- name:

  name of the manifold of interests, *"spd"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  Generate 5 observations; empirical covariance of normal observations.
#-------------------------------------------------------------------
#  Data Generation
d1 = array(0,c(3,3,5))
d2 = list()
for (i in 1:5){
  dat = matrix(rnorm(10*3),ncol=3)
  d1[,,i] = stats::cov(dat)
  d2[[i]] = d1[,,i]
}

#  Run
test1 = wrap.spd(d1)
test2 = wrap.spd(d2)
```
