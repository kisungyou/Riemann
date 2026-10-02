# Test of Uniformity on Sphere

Given \\N\\ observations \\\lbrace X_1, X_2, \ldots, X_M \brace\\ on
\\\mathcal{S}^{p-1}\\, it tests whether the data is distributed
uniformly on the sphere.

## Usage

``` r
sphere.utest(spobj, method = c("Rayleigh", "RayleighM"))
```

## Arguments

- spobj:

  a S3 `"riemdata"` class for \\N\\ Sphere-valued data.

- method:

  (case-insensitive) name of the test method containing

  `"Rayleigh"`

  :   original Rayleigh statistic.

  `"RayleighM"`

  :   modified Rayleigh statistic with better order of error.

## Value

a (list) object of `S3` class `htest` containing:

- statistic:

  a test statistic.

- p.value:

  \\p\\-value under \\H_0\\.

- alternative:

  alternative hypothesis.

- method:

  name of the test.

- data.name:

  name(s) of provided sample data.

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

Mardia KV, Jupp PE (eds.) (1999). *Directional Statistics*, Wiley Series
in Probability and Statistics. John Wiley and Sons, Inc., Hoboken, NJ,
USA. ISBN 978-0-470-31697-9 978-0-471-95333-3.

## See also

[`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)

## Examples

``` r
#-------------------------------------------------------------------
#   Compare Rayleigh's original and modified versions of the test
#-------------------------------------------------------------------
#  Data Generation
myobj = sphere.runif(n=100, p=5, type="riemdata")

#  Compare 2 versions : Original vs Modified Rayleigh
sphere.utest(myobj, method="rayleigh")
#> 
#>  Rayleigh Test of Uniformity on Sphere
#> 
#> data:  myobj
#> statistic = 4.8459, p-value = 0.435
#> alternative hypothesis: data is not uniformly distributed on 4-sphere.
#> 
sphere.utest(myobj, method="rayleighm")
#> 
#>  Modified Rayleigh Test of Uniformity on Sphere
#> 
#> data:  myobj
#> statistic = 4.8385, p-value = 0.4359
#> alternative hypothesis: data is not uniformly distributed on 4-sphere.
#> 

```
