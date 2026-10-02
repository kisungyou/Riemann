# Prepare Data on Grassmann Manifold

Grassmann manifold \\Gr(k,p)\\ is the set of \\k\\-planes, or
\\k\\-dimensional subspaces in \\R^p\\, which means that for a given
matrix \\Y \in \mathbf{R}{p\times k}\\, the column space \\SPAN(Y)\\ is
an element in Grassmann manifold. We use a convention that each element
in \\Gr(k,p)\\ is represented as an orthonormal basis (ONB) \\X \in
\mathbf{R}^{p\times k}\\ where \$\$X^\top X = I_k.\$\$ If not provided
in such a form, this wrapper takes a QR decomposition of the given data
to recover a corresponding ONB.

## Usage

``` r
wrap.grassmann(input)
```

## Arguments

- input:

  data matrices to be wrapped as `riemdata` class. Following inputs are
  considered,

  array

  :   an \\(p\times k\times n)\\ array where each slice along 3rd
      dimension is a \\k\\-subspace basis in dimension \\p\\.

  list

  :   a length-\\n\\ list whose elements are \\(p\times k)\\ basis for
      \\k\\-subspace.

## Value

a named `riemdata` S3 object containing

- data:

  a list of k-subspace basis matrices.

- size:

  size of each k-subspace basis matrix.

- name:

  name of the manifold of interests, *"grassmann"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  Generate 5 observations in Gr(2,4)
#-------------------------------------------------------------------
#  Generation
d1 = array(0,c(4,2,5))
d2 = list()
for (i in 1:5){
  d1[,,i] = matrix(rnorm(4*2), ncol=2)
  d2[[i]] = d1[,,i]
}

#  Run
test1 = wrap.grassmann(d1)
test2 = wrap.grassmann(d2)
```
