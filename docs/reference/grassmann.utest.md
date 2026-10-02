# Test of Uniformity on Grassmann Manifold

Given the data on Grassmann manifold \\Gr(k,p)\\, it tests whether the
data is distributed uniformly.

## Usage

``` r
grassmann.utest(grobj, method = c("Bing", "BingM"))
```

## Arguments

- grobj:

  a S3 `"riemdata"` class of Grassmann-valued data.

- method:

  (case-insensitive) name of the test method containing

  `"Bing"`

  :   Bingham statistic.

  `"BingM"`

  :   modified Bingham statistic with better order of error.

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

## Details

This is an asymptotic moment-based test, requiring independent
identically distributed observations. It need not detect every
nonuniform distribution. The modified expansion is unavailable for
ambient dimension two, where its displayed coefficients are singular;
full-dimensional subspaces are also excluded. Small-sample calibration
and the modified expansion remain experimental; a nominal chi-squared
p-value does not certify finite-sample size control.

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

[`wrap.grassmann`](https://www.kisungyou.com/Riemann/reference/wrap.grassmann.md)

## Examples

``` r
#-------------------------------------------------------------------
#   Compare Bingham's original and modified versions of the test
# 
# Test 1. sample uniformly from Gr(2,4)
# Test 2. use perturbed principal components from 'iris' data in R^4
#         which is concentrated around a point to reject H0.
#-------------------------------------------------------------------
## Data Generation
#  1. uniform data
myobj1 = grassmann.runif(n=100, k=2, p=4)

#  2. perturbed principal components
data(iris)
irdat = list()
for (n in 1:100){
   tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
   irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
}
myobj2 = wrap.grassmann(irdat)

## Test 1 : uniform data
grassmann.utest(myobj1, method="Bing")
#> 
#>  Bingham Test of Uniformity on Grassmann Manifold
#> 
#> data:  myobj1
#> statistic = 10.278, p-value = 0.3284
#> alternative hypothesis: data is not uniformly distributed on Gr(2,4).
#> 
grassmann.utest(myobj1, method="BingM")
#> 
#>  Modified Bingham Test of Uniformity on Grassmann Manifold
#> 
#> data:  myobj1
#> statistic = 10.275, p-value = 0.3287
#> alternative hypothesis: data is not uniformly distributed on Gr(2,4).
#> 

## Tests : iris data
grassmann.utest(myobj2, method="bINg")   # method names are 
#> 
#>  Bingham Test of Uniformity on Grassmann Manifold
#> 
#> data:  myobj2
#> statistic = 203.24, p-value < 2.2e-16
#> alternative hypothesis: data is not uniformly distributed on Gr(2,4).
#> 
grassmann.utest(myobj2, method="BiNgM")  # CASE - INSENSITIVE !
#> 
#>  Modified Bingham Test of Uniformity on Grassmann Manifold
#> 
#> data:  myobj2
#> statistic = 221, p-value < 2.2e-16
#> alternative hypothesis: data is not uniformly distributed on Gr(2,4).
#> 
```
