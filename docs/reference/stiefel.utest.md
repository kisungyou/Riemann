# Test of Uniformity on Stiefel Manifold

Given the data on Stiefel manifold \\St(k,p)\\, it tests whether the
data is distributed uniformly.

## Usage

``` r
stiefel.utest(stobj, method = c("Rayleigh", "RayleighM"))
```

## Arguments

- stobj:

  a S3 `"riemdata"` class for \\N\\ Stiefel-valued data.

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

[`wrap.stiefel`](https://www.kisungyou.com/Riemann/reference/wrap.stiefel.md)

## Examples

``` r
#-------------------------------------------------------------------
#   Compare Rayleigh's original and modified versions of the test
# 
# Test 1. sample uniformly from St(2,4)
# Test 2. use perturbed principal components from 'iris' data in R^4
#         which is concentrated around a point to reject H0.
#-------------------------------------------------------------------
## DATA GENERATION
#  1. uniform data
myobj1 = stiefel.runif(n=100, k=2, p=4)

#  2. perturbed principal components
data(iris)
irdat = list()
for (n in 1:100){
   tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
   irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
}
myobj2 = wrap.stiefel(irdat)

## TEST
#  1. uniform data
stiefel.utest(myobj1, method="Rayleigh")
#> 
#>  Rayleigh Test of Uniformity on Stiefel Manifold
#> 
#> data:  myobj1
#> statistic = 5.4403, p-value = 0.7096
#> alternative hypothesis: data is not uniformly distributed on St(2,4).
#> 
stiefel.utest(myobj1, method="RayleighM")
#> 
#>  Modified Rayleigh Test of Uniformity on Stiefel Manifold
#> 
#> data:  myobj1
#> statistic = 5.4279, p-value = 0.711
#> alternative hypothesis: data is not uniformly distributed on St(2,4).
#> 

#  2. concentrated data
stiefel.utest(myobj2, method="rayleIgh")   # method names are 
#> 
#>  Rayleigh Test of Uniformity on Stiefel Manifold
#> 
#> data:  myobj2
#> statistic = 52.309, p-value = 1.467e-08
#> alternative hypothesis: data is not uniformly distributed on St(2,4).
#> 
stiefel.utest(myobj2, method="raYleiGhM")  # CASE - INSENSITIVE !
#> 
#>  Modified Rayleigh Test of Uniformity on Stiefel Manifold
#> 
#> data:  myobj2
#> statistic = 53.415, p-value = 8.961e-09
#> alternative hypothesis: data is not uniformly distributed on St(2,4).
#> 
```
