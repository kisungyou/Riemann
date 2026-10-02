# Migrating analyses and saved models

The consolidation release changes numerical answers where the previous
implementation violated its mathematical contract. Keep the package
version and geometry with a saved analysis; retaining old fields does
not imply retaining an incorrect answer.

## Existing names and result fields

| Existing access | Current behavior |
|----|----|
| `riem.mean(x)$mean` | Preserved; convergence and objective fields are added. |
| `riem.pga(x)$embed` | Preserved; `ndim` is respected up to numerical rank. |
| `riem.kmeans(x)$cluster`, `$means`, `$score` | Preserved; starts and termination are recorded. |
| `predict(regression, newdata)` | Uses the fitted geometry automatically. |
| `geometry = "intrinsic"` on SPD | Resolves to `"affine_invariant"`. |
| `geometry = "extrinsic"` on SPD | Resolves to `"log_euclidean"`. |

The
[`riem.kmeans()`](https://www.kisungyou.com/Riemann/reference/riem.kmeans.md)
default is now Lloyd iteration. Set `algorithm = "MacQueen"` explicitly
if that sequential update policy is desired. Its manifold implementation
uses the shared Fréchet mean solver; approximate or nonconvex mean
solves do not receive a universal monotonicity guarantee.

## Numerical controls are honored

Invalid positive-integer controls now produce errors instead of silently
rounding or increasing the requested value. Short iteration budgets can
produce an unconverged estimate, which must be inspected before reuse.

``` r

x <- wrap.spd(list(diag(2), 4 * diag(2)))
mean_fit <- riem.mean(x, geometry = "affine_invariant", maxiter = 20, eps = 1e-9)
stopifnot(mean_fit$converged)
mean_fit$mean
#>      [,1] [,2]
#> [1,]    2    0
#> [2,]    0    2
stopifnot(max(abs(mean_fit$mean - 2 * diag(2))) < 1e-7)
invalid <- tryCatch(riem.mean(x, maxiter = 0), error = identity)
stopifnot(inherits(invalid, "error"))
conditionMessage(invalid)
#> [1] "maxiter must be a finite positive integer."
```

The mean of `I` and `4 I` is `2 I` under the affine-invariant metric.
Earlier iteration-scale errors could cycle in this case. The regression
suite now includes independent mathematical answers, together with
checks of actual implementation paths.

## Old regression objects without geometry

An old object may contain observations, responses, and a bandwidth but
lack the geometry used to fit it. These values cannot identify the
missing choice. Refit using the original recorded geometry. For a
one-time prediction, an explicit original geometry is accepted with a
warning; it is the caller’s responsibility to recover that information
from the analysis record.

``` r

x <- wrap.spd(lapply(c(1, 2, 4, 8), function(a) diag(c(a, 1))))
y <- c(0, 1, 0.5, 2)
fit <- riem.m2skreg(x, y, geometry = "log_euclidean", bandwidth = 0.5)
legacy <- fit
legacy$geometry <- NULL
legacy$schema_version <- NULL
missing_geometry <- tryCatch(predict(legacy, x), error = identity)
stopifnot(inherits(missing_geometry, "error"))
conditionMessage(missing_geometry)
#> [1] "This legacy fit has no geometry metadata. Supply its original geometry explicitly or refit it."
recovered <- predict(legacy, x, geometry = "log_euclidean")
#> Warning: Legacy fit: using the explicitly supplied geometry; the original
#> geometry cannot be verified.
stopifnot(isTRUE(all.equal(recovered, predict(fit, x))))
```

The warning in this example is deliberate. A saved modern model rejects
a different geometry at prediction time; refit explicitly to compare
another metric.

``` r

changed_geometry <- tryCatch(predict(fit, x, geometry = "affine_invariant"),
                              error = identity)
stopifnot(inherits(changed_geometry, "error"))
conditionMessage(changed_geometry)
#> [1] "Prediction geometry must match the fitted geometry; refit to change it."
```

## Refit old tangent and clustering objects

An old PGA result with only `center` and `embed` does not contain the
loadings, coordinate convention, or offset needed to transform new
observations. Refit it from the original training observations.
Similarly, old clustering lists should be refitted before using the new
prediction interface; their geometry cannot be reliably inferred from
centers alone.

``` r

training <- wrap.euclidean(rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 2)))
pca <- riem.pga(training, ndim = 2)
newdata <- wrap.euclidean(matrix(c(0.5, 0.8), 1))
predict(pca, newdata)
#>             PC1        PC2
#> [1,] 0.04832498 -0.0128334
riem.reconstruct(pca, predict(pca, newdata))$data[[1]]
#>      [,1]
#> [1,]  0.5
#> [2,]  0.8
```

Retain the order and labels of features, channels, and landmarks. The
wrapper and fitted contracts check recorded names, but they cannot
detect an undocumented semantic change to an unnamed column.

## Domain changes are explicit

Indefinite matrices are not SPD observations. Ill-conditioning and
positive definiteness are separate issues; adding a ridge or clipping
eigenvalues changes the data and belongs in the recorded preprocessing.
Sphere zero vectors cannot be normalized. Landmark analysis removes
location and scale and uses an orthogonal quotient that permits
reflections; verify that this matches the intended question.

Some additional manifold families remain available with an
`audit_pending` capability status. That label is not a claim of
comprehensive validation. Restrict scientific conclusions to the
operations and domains actually checked in the analysis, and retain
diagnostics and session information.

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Riemann_0.2.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1     hdrcde_3.5.0         dplyr_1.2.1         
#>  [4] farver_2.1.2         bitops_1.1-0         S7_0.2.2            
#>  [7] RCurl_1.98-1.20      fastmap_1.2.0        pracma_2.4.6        
#> [10] maotai_0.3.0         RANN_2.6.3           digest_0.6.39       
#> [13] lifecycle_1.0.5      cluster_2.1.8.3      rstiefel_1.0.1      
#> [16] magrittr_2.0.5       dbscan_1.2.6         compiler_4.5.2      
#> [19] rlang_1.3.0          sass_0.4.10          tools_4.5.2         
#> [22] mclustcomp_0.3.5     yaml_2.3.12          knitr_1.51          
#> [25] htmlwidgets_1.6.4    scatterplot3d_0.3-45 mclust_6.1.3        
#> [28] RColorBrewer_1.1-3   rainbow_3.8          KernSmooth_2.23-27  
#> [31] fda_6.3.0            Rtsne_0.17           desc_1.4.3          
#> [34] pcaPP_2.0-5          grid_4.5.2           clarabel_0.11.3     
#> [37] colorspace_2.1-3     ADMM_0.3.4           T4cluster_0.1.4     
#> [40] fastcluster_1.3.0    ggplot2_4.0.3        scales_1.4.0        
#> [43] iterators_1.0.14     MASS_7.3-66          mvtnorm_1.4-2       
#> [46] cli_3.6.6            rmarkdown_2.31       ragg_1.5.2          
#> [49] generics_0.1.4       otel_0.2.0           RSpectra_0.16-2     
#> [52] RcppDE_0.1.9         CVXR_1.9.2           scs_3.2.7           
#> [55] fds_1.9              cachem_1.1.0         splines_4.5.2       
#> [58] parallel_4.5.2       vctrs_0.7.3          T4transport_0.1.9   
#> [61] Matrix_1.7-6         jsonlite_2.0.0       systemfonts_1.3.2   
#> [64] foreach_1.5.2        gsignal_0.3-7        jquerylib_0.1.4     
#> [67] glue_1.8.1           pkgdown_2.2.1        codetools_0.2-20    
#> [70] DEoptim_2.2-8        gtable_0.3.6         osqp_1.0.0          
#> [73] gmp_0.7-5.1          tibble_3.3.1         pillar_1.11.1       
#> [76] htmltools_0.5.9      deSolve_1.42         R6_2.6.1            
#> [79] textshaping_1.0.5    Rdpack_2.6.6         ks_1.15.3           
#> [82] doParallel_1.0.17    lpSolve_5.6.23       evaluate_1.0.5      
#> [85] lattice_0.23-1       rbibutils_2.4.1      backports_1.5.1     
#> [88] highs_1.14.0-2       bslib_0.12.0         Rcpp_1.1.2          
#> [91] checkmate_2.3.4      xfun_0.60            fs_2.1.0            
#> [94] Rdimtools_1.1.5      pkgconfig_2.0.3
```
