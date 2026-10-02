# Geometry-consistent statistical workflows

A geometry specifies the distances and operations used by an analysis.
Select it before fitting, inspect the numerical result, and reuse the
fitted object for new observations. The examples below use small
synthetic observations. They illustrate software behavior, without
making scientific claims about an application dataset.

## Validate observations and select a geometry

Construct positive-definite matrices by exponentiating real symmetric
matrices. The off-diagonal chart entries give noncommuting observations.

``` r

matrix_exp <- function(a) {
  e <- eigen(a, symmetric = TRUE)
  e$vectors %*% diag(exp(e$values)) %*% t(e$vectors)
}
t <- seq(-0.6, 0.6, length.out = 12)
matrices <- lapply(t, function(a) {
  matrix_exp(matrix(c(a, 0.2 * sin(3 * a), 0.2 * sin(3 * a), -0.5 * a), 2))
})
x <- wrap.spd(matrices)
airm <- riem.geometry(x, "affine_invariant")
log_euclidean <- riem.geometry(x, "log_euclidean")
airm
#> Geometry: affine_invariant on spd 
#> Representation: 2 x 2 ; intrinsic dimension: 3 
#> Validation: core_contract
log_euclidean
#> Geometry: log_euclidean on spd 
#> Representation: 2 x 2 ; intrinsic dimension: 3 
#> Validation: core_contract
```

Affine-invariant and log-Euclidean distances describe different
problems. The legacy SPD names `"intrinsic"` and `"extrinsic"` resolve
to these two choices, respectively. The capability table records which
operation routes are available; its validation status also identifies
families still awaiting a broader audit.

``` r

capability <- riem.capabilities()
capability[capability$manifold_id %in% c("spd", "sphere", "landmark"), ]
#>    manifold_id   backend        geometry_id distance  mean median tangent
#> 1          spd intrinsic   affine_invariant     TRUE  TRUE   TRUE    TRUE
#> 2       sphere intrinsic              round     TRUE  TRUE   TRUE    TRUE
#> 4     landmark intrinsic   shape_orthogonal     TRUE  TRUE   TRUE    TRUE
#> 11         spd extrinsic      log_euclidean     TRUE  TRUE   TRUE    TRUE
#> 12      sphere extrinsic            chordal     TRUE  TRUE   TRUE   FALSE
#> 14    landmark extrinsic procrustes_chordal     TRUE FALSE  FALSE   FALSE
#>    validation_status          median_estimand
#> 1      core_contract           frechet_median
#> 2      core_contract           frechet_median
#> 4      audit_pending           frechet_median
#> 11     core_contract           frechet_median
#> 12     core_contract projected_ambient_median
#> 14     audit_pending projected_ambient_median
```

Positive definiteness is checked independently of full rank. This
deliberately invalid example must produce an error; the vignette
verifies that the error occurs.

``` r

invalid <- tryCatch(wrap.spd(list(diag(c(-1, 1)))), error = identity)
stopifnot(inherits(invalid, "error"))
conditionMessage(invalid)
#> [1] "wrap.spd: observation 1 must be a finite real symmetric positive-definite matrix."
```

## Inspect a numerical summary

The shared solver records its objective and termination state. A finite
matrix by itself is not a convergence certificate.

``` r

m <- riem.mean(x, geometry = airm, maxiter = 200, eps = 1e-8, trace = TRUE)
m
#> Geometric mean under affine_invariant on spd 
#> Termination: stationary ; converged: TRUE ; iterations: 4 
#> Objective: 0.224607 
#> Estimand: intrinsic_frechet_mean 
#>               [,1]          [,2]
#> [1,]  1.000000e+00 -1.663168e-10
#> [2,] -1.663168e-10  1.000000e+00
stopifnot(m$converged)
c(objective = m$objective, iterations = m$iterations)
#>  objective iterations 
#>   0.224607   4.000000
mlog <- riem.mean(x, geometry = log_euclidean)
stopifnot(mlog$converged)
riem.median(x, geometry = log_euclidean, maxiter = 200)
#> Geometric median under log_euclidean on spd 
#> Termination: stationary ; converged: TRUE ; iterations: 0 
#> Objective: 0.41860924 
#> Estimand: extrinsic_frechet_median 
#>               [,1]          [,2]
#> [1,]  1.000000e+00 -2.489623e-18
#> [2,] -2.489623e-18  1.000000e+00
```

## Fit, predict, and reconstruct tangent components

Tangent PCA approximates variation near a trained reference. It retains
a metric-preserving coordinate convention and any coordinate offset. Its
explained variance is a property of this tangent representation, rather
than exact global variance on the manifold.

``` r

pca <- riem.pga(x, ndim = 2, geometry = log_euclidean)
pca
#> Tangent PCA on spd ( log_euclidean )
#> 12 observations; 2 components retained; numerical rank 2
summary(pca)
#> Variance explained in the fitted tangent representation 
#>                  PC1         PC2
#> variance   0.2435162 0.001509638
#> proportion 0.9938389 0.006161138
#> cumulative 0.9938389 1.000000000
stopifnot(isTRUE(all.equal(unname(predict(pca, x)), unname(pca$embed))))
new_x <- wrap.spd(list(matrix_exp(matrix(c(0.15, 0.08, 0.08, -0.1), 2))))
new_scores <- predict(pca, new_x)
new_scores
#>            PC1        PC2
#> [1,] 0.2107759 0.01932642
approximation <- riem.reconstruct(pca, new_scores)
approximation$data[[1]]
#>            [,1]       [,2]
#> [1,] 1.17698470 0.08355381
#> [2,] 0.08355381 0.92632326
stopifnot(min(eigen(approximation$data[[1]], symmetric = TRUE)$values) > 0)
```

Prediction and reconstruction use the saved reference. With fewer
components than the data’s rank, reconstruction is an approximation.
Requests above numerical rank are truncated with a warning and the
effective dimension is recorded.

The sphere uses an orthonormal tangent basis. The following observations
lie in a small neighborhood where the reference and logarithms are well
defined.

``` r

sphere <- wrap.sphere(cbind(rep(1, 10), matrix(rnorm(20, sd = 0.08), 10, 2)))
sphere_pca <- riem.pga(sphere, ndim = 2, geometry = "round")
reconstructed <- riem.reconstruct(sphere_pca)
stopifnot(all(vapply(reconstructed$data, function(a) abs(sum(a^2) - 1) < 1e-10,
                    logical(1))))
plot(sphere_pca)
```

![](geometry-workflows_files/figure-html/sphere-pca-1.png)

Antipodal logarithms are nonunique. Unsupported geometry-operation
combinations also produce explicit errors rather than silently switching
metrics.

``` r

unsupported <- tryCatch(riem.pga(sphere, geometry = "chordal"), error = identity)
stopifnot(inherits(unsupported, "error"))
conditionMessage(unsupported)
#> [1] "Geometry 'chordal' does not support tangent with a compatible implementation."
```

## Reuse a clustering model

The default Lloyd driver uses the shared mean solver and records each
start. Labels returned by prediction refer to the final trained centers.

``` r

set.seed(1702)
clusters <- riem.kmeans(x, k = 2, geometry = log_euclidean, nstart = 2)
clusters
#> Manifold k-means: 2 clusters on log_euclidean 
#> Within-cluster sum of squares: 0.6083241 ; termination: stable_assignment
clusters$starts
#>   start valid converged       termination objective iterations empty_clusters
#> 1     1  TRUE      TRUE stable_assignment 0.6083241          3              0
#> 2     2  TRUE      TRUE stable_assignment 0.6083241          1              0
#>   error
#> 1      
#> 2
stopifnot(identical(predict(clusters, x), clusters$cluster))
predict(clusters, new_x)
#> [1] 1
```

## Tune scalar-response regression using every held-out fold

Here the predictor is a matrix and the response is numeric. Geometry is
fixed before tuning. The fold assignments are stored and every
observation contributes once to the held-out sum of squared errors. For
grouped or temporal data, supply splits appropriate to the sampling
design. If preprocessing is learned from data, it must be repeated
within the corresponding training folds.

``` r

y <- sin(2 * t) + rnorm(length(t), sd = 0.03)
regression <- riem.m2skregCV(
  x, y, geometry = log_euclidean, bandwidths = c(0.1, 0.25, 0.5),
  foldid = rep(1:3, length.out = length(t))
)
regression$errors
#>      bandwidth        SSE
#> [1,]      0.10 0.05392788
#> [2,]      0.25 0.14691578
#> [3,]      0.50 0.99031984
regression$fold_errors
#>               1           2          3
#> [1,] 0.01809862 0.007890027 0.02793924
#> [2,] 0.03641123 0.029304317 0.08120023
#> [3,] 0.29816855 0.264589651 0.42756163
prediction <- predict(regression, new_x, diagnostics = TRUE)
prediction
#> $prediction
#> [1] 0.2900565
#> 
#> $diagnostics
#>   effective_n nearest_distance nearest_over_bandwidth
#> 1    2.408551       0.03040861              0.3040861
stopifnot(identical(regression$geometry$geometry_id, "log_euclidean"))
```

A finite prediction can still have weak local support. Effective
neighborhood size and distance relative to bandwidth help assess that
limitation; numerical stability does not establish statistical support.

## Preserve the fitted state

``` r

path <- tempfile(fileext = ".rds")
saveRDS(pca, path)
restored <- readRDS(path)
unlink(path)
stopifnot(isTRUE(all.equal(predict(restored, new_x), new_scores)))
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

Regular landmark tangent workflows use equivalence under the full
orthogonal group, including reflections, and align new observations to
the trained reference. Singular alignments and reconstructions outside
the supported local chart are rejected. This convention must agree with
the scientific meaning of the landmarks; it is not an
orientation-preserving shape convention.
