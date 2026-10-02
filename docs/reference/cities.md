# Data : Populated Cities in the U.S.

As of January 2006, there are 60 cities in the contiguous U.S. with
population size larger than \\300000\\. We extracted information of the
cities from the data delivered by maps package. Variables `coord` and
`cartesian` are two identical representations of locations, which can be
mutually converted by
[`sphere.convert`](https://www.kisungyou.com/Riemann/reference/sphere.convert.md).

## Usage

``` r
data(cities)
```

## Format

a named list containing

- names:

  a length-\\60\\ vector of city names.

- coord:

  a \\(60\times 2)\\ matrix of latitude and longitude.

- cartesian:

  a \\(60\times 3)\\ matrix of cartesian coordinates on the unit sphere.

- population:

  a length-\\60\\ vector of cities' populations.

## See also

[`wrap.sphere`](https://www.kisungyou.com/Riemann/reference/wrap.sphere.md)

## Examples

``` r
# \donttest{
## LOAD THE DATA AND WRAP AS RIEMOBJ
data(cities)
myriem = wrap.sphere(cities$cartesian)

## COMPUTE INTRINSIC/EXTRINSIC MEANS
intmean = as.vector(riem.mean(myriem, geometry="intrinsic")$mean)
extmean = as.vector(riem.mean(myriem, geometry="extrinsic")$mean)

## CONVERT TO GEOGRAPHIC COORDINATES (Lat/Lon)
geo.int = sphere.xyz2geo(intmean)
geo.ext = sphere.xyz2geo(extmean)
# }
```
