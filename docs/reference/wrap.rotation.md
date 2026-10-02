# Prepare Data on Rotation Group

Rotation group, also known as special orthogonal group, is a Riemannian
manifold \$\$SO(p) = \lbrace Q \in \mathbf{R}^{p\times p}~\vert~ Q^\top
Q = I, \textrm{det}(Q)=1 \rbrace \$\$ where the name originates from an
observation that when \\p=2,3\\ these matrices are rotation of
shapes/configurations.

## Usage

``` r
wrap.rotation(input)
```

## Arguments

- input:

  data matrices to be wrapped as `riemdata` class. Following inputs are
  considered,

  array

  :   a \\(p\times p\times n)\\ array where each slice along 3rd
      dimension is a rotation matrix.

  list

  :   a length-\\n\\ list whose elements are \\(p\times p)\\ rotation
      matrices.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times p)\\ rotation matrices.

- size:

  size of each rotation matrix.

- name:

  name of the manifold of interests, *"rotation"*

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#-------------------------------------------------------------------
## DATA GENERATION
d1 = array(0,c(3,3,5))
d2 = list()
for (i in 1:5){
  single  = qr.Q(qr(matrix(rnorm(9),nrow=3)))
  d1[,,i] = single
  d2[[i]] = single
}

## RUN
test1 = wrap.rotation(d1)
test2 = wrap.rotation(d2)
```
