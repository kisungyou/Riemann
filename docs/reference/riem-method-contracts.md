# Method Contracts and Validation Status

Riemann distinguishes validated input representations, available
geometry primitives, algorithm contracts, and statistical calibration. A
wrapped point or available distance does not establish every algorithm's
assumptions.

## Details

The installed file `method-contracts.csv` contains one disposition for
every exported function and registered S3 method, including its domain,
assumptions, restrictions, and evidence scope. Read it with the example
below. The four possible dispositions are:

- core:

  The consolidated input, geometry, summary, clustering, tangent PCA, or
  scalar regression contract, within its documented capabilities.

- validated_restricted:

  A specific legacy contract checked against independent formulas or
  reference fixtures; the stated restrictions apply.

- experimental:

  A retained legacy interface whose full numerical, mathematical, or
  statistical contract has not been independently verified.

- disabled:

  An unavailable interface with a documented incompatible
  implementation. Current geometry restrictions can disable combinations
  without disabling an entire exported function.

A test-file reference is an evidence pointer, not branch coverage or a
proof. No disposition guarantees a global optimum, an identifiable
model, finite sample calibration, or validity for every manifold.
Experimental methods are retained for compatibility and exploration;
their source references and examples should not be read as completed
validation.

Common distinctions are particularly relevant: nonnegative graph
affinities or regression weights need not form a positive-semidefinite
kernel; a metric dissimilarity need not be Euclidean; a local tangent
representation need not solve exact nonlinear principal geodesic
analysis; and a permutation test requires exchangeability under its
specified null. Rayleigh/Bingham and Frechet mean/variance tests need
not detect all distributional alternatives.

[`riem.capabilities()`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md)
records primitive availability and geometry status. The table here
separately records method status. The package maintenance inventory
records source and dependency ownership; the release/check logs supply
actual execution evidence for each environment.

## See also

[`riem.capabilities`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md),
[`riem.geometry`](https://www.kisungyou.com/Riemann/reference/riem.geometry.md)

## Examples

``` r
contracts <- read.csv(system.file("method-contracts.csv", package = "Riemann"))
contracts[contracts$name %in% c("riem.kpca", "riem.fanova"),
          c("name", "disposition", "limitations")]
#>           name          disposition
#> 44 riem.fanova validated_restricted
#> 56   riem.kpca validated_restricted
#>                                                                                                                                                                                       limitations
#> 44             Not an omnibus all-distribution test; paired/clustered data require another resampling design. Empirical nondegeneracy and solver convergence do not prove population assumptions.
#> 56 Empirical negative spectrum is reported, warns and is explicitly truncated or rejected. Truncation is an approximation; an empirically PSD sample does not prove global positive definiteness.
```
