# Supported Geometries on SPD Manifold

SPD manifold is a well-studied space in that there have been many
geometries proposed on the space. For special functions on under SPD
category, this function finds whether there exists a matching name that
is currently supported in Riemann. If there is none, it will return an
error message.

## Usage

``` r
spd.geometry(geometry)
```

## Arguments

- geometry:

  name of supported geometries, including

  AIRM

  :   Affine-Invariant Riemannian Metric.

  LERM

  :   Log-Euclidean Riemannian Metric.

  Jeffrey

  :   Jeffrey's divergence.

  Stein

  :   Stein's metric.

  Wasserstein

  :   2-Wasserstein geometry.

## Value

a matching name in lower-case.

## Validation status

This retained legacy interface is experimental. Its full numerical and
statistical contract has not been independently verified across
supported inputs. See
[`riem-method-contracts`](https://www.kisungyou.com/Riemann/reference/riem-method-contracts.md)
and the installed contract table for method-specific assumptions,
restrictions, and evidence scope.

## Examples

``` r
# it just returns a small-letter string.
mygeom = spd.geometry("stein")
```
