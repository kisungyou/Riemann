# Generate Uniform Samples on Sphere

It generates \\n\\ random samples from \\\mathcal{S}^{p-1}\\. For
convenient usage of users, we provide a number of options in terms of
the return type.

## Usage

``` r
sphere.runif(n, p, type = c("list", "matrix", "riemdata"))
```

## Arguments

- n:

  number of samples to be generated.

- p:

  original dimension (of the ambient space).

- type:

  return type;

  `"list"`

  :   a length-\\n\\ list of length-\\p\\ vectors.

  `"matrix"`

  :   a \\(n\times p)\\ where rows are unit vectors.

  `"riemdata"`

  :   a S3 object. See
      [`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)
      for more details (*Default*).

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

[`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)

## Examples

``` r
#-------------------------------------------------------------------
#                       Draw Samples on Sphere
#
# Multiple return types on S^4 in R^5
#-------------------------------------------------------------------
dat.list = sphere.runif(n=10, p=5, type="list")
dat.matx = sphere.runif(n=10, p=5, type="matrix")
dat.riem = sphere.runif(n=10, p=5, type="riemdata")
```
