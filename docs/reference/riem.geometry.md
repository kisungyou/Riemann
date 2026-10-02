# Inspect and Resolve a Statistical Geometry

Geometry is resolved once when a model is fitted and retained for
prediction. For SPD matrices, `"intrinsic"` means `"affine_invariant"`
and `"extrinsic"` means `"log_euclidean"`. On the sphere the
corresponding names are `"round"` and `"chordal"`. Euclidean aliases
coincide. Landmark geometry currently identifies configurations under
the full orthogonal group, including reflections. A capability records
an available operation; its validation status must also be considered
before interpreting results.

## Usage

``` r
riem.geometry(riemobj, geometry = NULL)

riem.capabilities()
```

## Arguments

- riemobj:

  A wrapped `riemdata` object.

- geometry:

  A geometry name, a saved `riem_geometry` specification, or `NULL` for
  the default intrinsic geometry.

## Value

`riem.geometry` returns a serializable geometry specification.
`riem.capabilities` returns the operation registry as a data frame.

## Examples

``` r
x <- wrap.spd(list(diag(2), 2 * diag(2)))
riem.geometry(x, "log_euclidean")
#> Geometry: log_euclidean on spd 
#> Representation: 2 x 2 ; intrinsic dimension: 3 
#> Validation: core_contract 
riem.capabilities()
#>    manifold_id   backend                       geometry_id distance  mean
#> 1          spd intrinsic                  affine_invariant     TRUE  TRUE
#> 2       sphere intrinsic                             round     TRUE  TRUE
#> 3    euclidean intrinsic                         euclidean     TRUE  TRUE
#> 4     landmark intrinsic                  shape_orthogonal     TRUE  TRUE
#> 5    grassmann intrinsic                  principal_angles     TRUE  TRUE
#> 6      stiefel intrinsic     stiefel_intrinsic_unavailable    FALSE FALSE
#> 7     rotation intrinsic                rotation_frobenius     TRUE  TRUE
#> 8  multinomial intrinsic                        fisher_rao     TRUE  TRUE
#> 9         spdk intrinsic                 factor_procrustes     TRUE  TRUE
#> 10 correlation intrinsic  correlation_quotient_unavailable    FALSE FALSE
#> 11         spd extrinsic                     log_euclidean     TRUE  TRUE
#> 12      sphere extrinsic                           chordal     TRUE  TRUE
#> 13   euclidean extrinsic                         euclidean     TRUE  TRUE
#> 14    landmark extrinsic                procrustes_chordal     TRUE FALSE
#> 15   grassmann extrinsic                 projector_chordal     TRUE  TRUE
#> 16     stiefel extrinsic                     frame_chordal     TRUE  TRUE
#> 17    rotation extrinsic                  rotation_chordal     TRUE  TRUE
#> 18 multinomial extrinsic                      sqrt_chordal     TRUE  TRUE
#> 19        spdk extrinsic      factor_extrinsic_unavailable    FALSE FALSE
#> 20 correlation extrinsic correlation_extrinsic_unavailable    FALSE FALSE
#>    median tangent                validation_status          median_estimand
#> 1    TRUE    TRUE                    core_contract           frechet_median
#> 2    TRUE    TRUE                    core_contract           frechet_median
#> 3    TRUE    TRUE                    core_contract           frechet_median
#> 4    TRUE    TRUE                    audit_pending           frechet_median
#> 5    TRUE   FALSE                    audit_pending           frechet_median
#> 6   FALSE   FALSE restricted_inconsistent_geometry           frechet_median
#> 7    TRUE   FALSE                    audit_pending           frechet_median
#> 8    TRUE   FALSE                    audit_pending           frechet_median
#> 9    TRUE   FALSE                    audit_pending           frechet_median
#> 10  FALSE   FALSE restricted_inconsistent_geometry           frechet_median
#> 11   TRUE    TRUE                    core_contract           frechet_median
#> 12   TRUE   FALSE                    core_contract projected_ambient_median
#> 13   TRUE    TRUE                    core_contract           frechet_median
#> 14  FALSE   FALSE                    audit_pending projected_ambient_median
#> 15   TRUE   FALSE                    audit_pending projected_ambient_median
#> 16   TRUE   FALSE                    audit_pending projected_ambient_median
#> 17   TRUE   FALSE                    audit_pending projected_ambient_median
#> 18   TRUE   FALSE                    audit_pending projected_ambient_median
#> 19  FALSE   FALSE restricted_inconsistent_geometry projected_ambient_median
#> 20  FALSE   FALSE restricted_inconsistent_geometry projected_ambient_median
```
