# Convert between Cartesian Coordinates and Geographic Coordinates

In geospatial data analysis, it is common to consider locations on the
Earth as data. These locations, usually provided by latitude and
longitude, are not directly applicable for spherical data analysis. We
provide two functions - `sphere.geo2xyz` and `sphere.xyz2geo` - that
convert geographic coordinates in longitude/latitude into a unit-norm
vector on \\\mathcal{S}^2\\, and vice versa. As a convention, latitude
and longitude are represented as *decimal degrees*.

## Usage

``` r
sphere.geo2xyz(lat, lon)

sphere.xyz2geo(xyz)
```

## Arguments

- lat:

  latitude (in decimal degrees).

- lon:

  longitude (in decimal degrees).

- xyz:

  a unit-norm vector in \\\mathcal{S}^{2}\\.

## Value

transformed data.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## Examples

``` r
## EXAMPLE DATA WITH POPULATED US CITIES
data(cities)

## SELECT ALBUQUERQUE
geo = cities$coord[1,]
xyz = cities$cartesian[1,]

## CHECK TWO INPUT TYPES AND THEIR CONVERSIONS
sphere.geo2xyz(geo[1], geo[2])
#> [1] -0.2339521 -0.7837773  0.5752908
sphere.xyz2geo(xyz)
#> [1]   35.12 -106.62
```
