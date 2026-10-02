# Prepare Data on Correlation Manifold

The collection of correlation matrices is considered as a subset (and
quotient) of the well-known SPD manifold. In our package, it is defined
as \$\$\mathcal{C}\_{++}^p = \lbrace X \in \mathbf{R}^{p\times p}
~\vert~ X^\top = X,~ X \succ 0,~ \textrm{diag}(X) = 1 \rbrace\$\$ where
\\X \succ 0\\ means strictly positive definite. Please note that the
geometry involving semi-definite correlation matrices is not the
objective here.

## Usage

``` r
wrap.correlation(input)
```

## Arguments

- input:

  correlation data matrices to be wrapped as `riemdata` class. Following
  inputs are considered,

  array

  :   an \\(p\times p\times n)\\ array where each slice along 3rd
      dimension is a correlation matrix.

  list

  :   a length-\\n\\ list whose elements are \\(p\times p)\\ correlation
      matrices.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times p)\\ correlation matrices.

- size:

  size of each correlation matrix.

- name:

  name of the manifold of interests, *"correlation"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  5 observations; empirical correlation of normal observations.
#-------------------------------------------------------------------
#  Data Generation
d1 = array(0,c(3,3,5))
d2 = list()
for (i in 1:5){
  dat = matrix(rnorm(10*3),ncol=3)
  d1[,,i] = stats::cor(dat)
  d2[[i]] = d1[,,i]
}

#  Run
test1 = wrap.correlation(d1)
test2 = wrap.correlation(d2)
```
