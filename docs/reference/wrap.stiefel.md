# Prepare Data on (Compact) Stiefel Manifold

Stiefel manifold \\St(k,p)\\ is the set of \\k\\-frames in
\\\mathbf{R}^p\\, which is indeed a Riemannian manifold. For usage in
Riemann package, each data point is represented as a matrix by the
convention \$\$St(k,p) = \lbrace X \in \mathbf{R}^{p\times k} ~\vert~
X^\top X = I_k \rbrace\$\$ which means that columns are orthonormal.
When the provided matrix is not an orthonormal basis as above,
`wrap.stiefel` applies orthogonalization to extract valid basis
information.

## Usage

``` r
wrap.stiefel(input)
```

## Arguments

- input:

  data matrices to be wrapped as `riemdata` class. Following inputs are
  considered,

  array

  :   a \\(p\times k\times n)\\ array where each slice along 3rd
      dimension is a \\k\\-frame.

  list

  :   a length-\\n\\ list whose elements are \\(p\times k)\\
      \\k\\-frames.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\k\\-frame orthonormal matrices.

- size:

  size of each \\k\\-frame basis matrix.

- name:

  name of the manifold of interests, *"stiefel"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#
#  Generate 5 observations in St(2,4)
#-------------------------------------------------------------------
#  Data Generation by QR Decomposition
d1 = array(0,c(4,2,5))
d2 = list()
for (i in 1:5){
  d1[,,i] = qr.Q(qr(matrix(rnorm(4*2),ncol=2)))
  d2[[i]] = d1[,,i]
}

#  Run
test1 = wrap.stiefel(d1)
test2 = wrap.stiefel(d2)
```
