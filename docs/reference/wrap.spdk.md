# Prepare Data on Positive Semidefinite Manifold of Fixed Rank

When \\(p\times p)\\ positive semidefinite matrices are of fixed rank
\\k \< p\\, they form a geometric structure represented by \\(p\times
k)\\ matrices, \$\$SPD(k,p) = \lbrace X \in \mathbf{R}^{(p\times
p)}~\vert~ Y Y^\top = X, \textrm{rank}(X) = k \rbrace\$\$ It's key
difference from \\\mathcal{S}\_{++}^p\\ is that all matrices should be
of fixed rank \\k\\ where \\k\\ is usually smaller than \\p\\. Inputs
are given as \\(p\times p)\\ matrices with specified \\k\\ and
`wrap.spdk` decomposes matrices whose numerical rank is exactly \\k\\.
Higher-rank inputs are rejected; rank reduction must be performed
explicitly before wrapping. Numerical rank uses eigenvalues exceeding
\\64p\\ times machine precision relative to the largest absolute
eigenvalue. Negative eigenvalues beyond this roundoff tolerance are
rejected.

## Usage

``` r
wrap.spdk(input, k)
```

## Arguments

- input:

  data matrices to be wrapped as `riemdata` class. Following inputs are
  considered,

  array

  :   a \\(p\times p\times n)\\ array where each slice along 3rd
      dimension is a rank-\\k\\ matrix.

  list

  :   a length-\\n\\ list whose elements are \\(p\times p)\\ matrices of
      rank-\\k\\.

- k:

  rank of the positive semidefinite matrices, an integer in
  \\\[1,p-1\]\\.

## Value

a named `riemdata` S3 object containing

- data:

  a list of \\(p\times k)\\ representation of the corresponding
  rank-\\k\\ SPSD matrices.

- size:

  size of each representation matrix.

- name:

  name of the manifold of interests, *"spdk"*

## References

Journée M, Bach F, Absil P, Sepulchre R (2010). “Low-rank optimization
on the cone of positive semidefinite matrices.” *SIAM Journal on
Optimization*, **20**(5), 2327–2351.

## Examples

``` r
#-------------------------------------------------------------------
#                 Checker for Two Types of Inputs
#-------------------------------------------------------------------
#  Data Generation
d1 = array(0,c(10,10,3))
d2 = list()
for (i in 1:3){
  dat = matrix(rnorm(10*2),nrow=10)
  d1[,,i] = tcrossprod(dat)
  d2[[i]] = d1[,,i]
}

#  Run
test1 = wrap.spdk(d1, k=2)
test2 = wrap.spdk(d2, k=2)

```
