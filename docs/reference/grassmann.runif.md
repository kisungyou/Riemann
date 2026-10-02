# Generate Uniform Samples on Grassmann Manifold

It generates \\n\\ random samples from Grassmann manifold \\Gr(k,p)\\.

## Usage

``` r
grassmann.runif(n, k, p, type = c("list", "array", "riemdata"))
```

## Arguments

- n:

  number of samples to be generated.

- k:

  dimension of the subspace.

- p:

  original dimension (of the ambient space).

- type:

  return type;

  `"list"`

  :   a length-\\n\\ list of \\(p\times k)\\ basis of \\k\\-subspaces.

  `"array"`

  :   a \\(p\times k\times n)\\ 3D array whose slices are \\k\\-subspace
      basis.

  `"riemdata"`

  :   a S3 object. See
      [`wrap.grassmann`](https://www.kisungyou.com/Riemann/reference/wrap.grassmann.md)
      for more details.

## Value

an object from one of the above by `type` option.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## References

Chikuse Y (2003). *Statistics on Special Manifolds*, volume 174 of
*Lecture Notes in Statistics*. Springer New York, New York, NY. ISBN
978-0-387-00160-9 978-0-387-21540-2.

## See also

[`stiefel.runif`](https://www.kisungyou.com/Riemann/reference/stiefel.runif.md),
[`wrap.grassmann`](https://www.kisungyou.com/Riemann/reference/wrap.grassmann.md)

## Examples

``` r
#-------------------------------------------------------------------
#                 Draw Samples on Grassmann Manifold 
#-------------------------------------------------------------------
#  Multiple Return Types with 3 Observations of 5-dim subspaces in R^10
dat.list = grassmann.runif(n=3, k=5, p=10, type="list")
dat.arr3 = grassmann.runif(n=3, k=5, p=10, type="array")
dat.riem = grassmann.runif(n=3, k=5, p=10, type="riemdata")
```
