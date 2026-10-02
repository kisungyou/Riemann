# Riemann 0.2.0

* Correct the October 2026 audit findings: Sammon curvature and descent,
  coreset sampling with replacement and weighted clustering, stable Stein and
  Bures-Wasserstein distances, spherical likelihood bounds and log densities,
  Stiefel QR signs, Grassmann optimizer returns, nearest-neighbor self exclusion,
  Grassmann metric-learning dimensions, and extreme-weight tangent variances.
  Coreset clustering now honors `maxiter`; coincident samples have an explicit
  zero-dispersion sampling rule. Help clarifies that first-order stationarity
  does not certify a local minimum. Independent regression fixtures accompany
  these corrections; legacy experimental status remains where appropriate.

* The September 21, 2026 source revision licenses author-owned package software
  under GPL-3 and retains third-party attribution. Earlier MIT releases and
  preserved archives keep their original licenses. Package source comments
  and roxygen text are now ASCII; formatted help retains accented names.

* The September 21 r3 documentation revision identifies the ERP matrices as
  xDAWN-augmented MEG covariances through a numerical reconstruction of the
  historical recipe. It records the MNE source's CC0 declaration and retains
  GPL-2 terms for the shapes-derived gorilla data. Data and algorithms are unchanged.

* Correct SPD validation and the affine-invariant metric. Symmetric matrix
  functions replace general complex-valued matrix logarithms in core SPD paths.
* Preserve tiny sphere angles and report nonunique antipodal logarithms and
  degenerate projected summaries explicitly.
* Share mean and median solvers between public summaries and internal callers,
  with validated weights/controls, descent, and actual termination diagnostics.
* Preserve fitted regression geometry, stabilize kernel weights, and correct
  all-fold cross-validation, singleton subsets, and candidate error tables.
* Add explicit geometry specifications and a capability registry. Existing
  intrinsic/extrinsic aliases remain available with their documented meanings.
* Tangent PCA now honors requested rank, preserves the chosen metric, and
  supports prediction and reconstruction under a trained reference.
* K-means uses validated Fréchet-center updates, per-start diagnostics, and
  fixed-center prediction. The default is now Lloyd; MacQueen remains available.
* Landmark wrapping explicitly uses O(p), including reflections, and rejects
  singular configurations. The default alignment reference is now the first
  normalized configuration without the former extra covariance-eigenvector
  rotation; quotient distances retain their interpretation.
* Correct additional Grassmann, rotation, positive-simplex and fixed-rank PSD
  primitives, with explicit cut-locus and projection-uniqueness checks. Intrinsic
  Stiefel and correlation operations are restricted pending consistent geometry.
* Reject invalid correlation matrices, negative simplex masses, rank-deficient
  frames and implicit fixed-rank truncation. Preserve singleton dimensions and
  named feature order through fitted summaries and prediction.
* Correct FANOVA scaling, permutation bookkeeping, kernel-PCA scores, classical
  scaling spectra and the mutual-neighbor Isomap graph. Updated help states the
  narrower assumptions and experimental status of remaining legacy routines.
* Cache log-Euclidean matrix transforms within pairwise and cross-distance calls.
  Add executed data workflows, matched software comparisons, and migration help.
* Validation and remaining publication requirements are tracked in the separate
  write-JSS-Riemann/maintenance/ workspace. Journal submission remains a separate step.

# Riemann 0.1.7

* Removed dependence on `CVXR` per recent API changes.

# Riemann 0.1.5

* Code covarage and unit testing start to be integrated.

# Riemann 0.1.4

* Two mixture models on the unit hypersphere added.
* Supports for spherical Laplace distribution added.
* `riem.phate()` is added for sub-manifold learning/visualization.
* New data from EEG ERPs is available by `data(ERP)`.
* `spd.wassbary()` is added to compute Wasserstein barycenter of Gaussian distributions.

# Riemann 0.1.3

* Modified `mle.spnorm()` for controls with user-defined stopping criteria.
* Fixed log-likelihood evaluation for `mixspnorm()`.

# Riemann 0.1.2

* Added `riem.m2skreg()` for manifold-to-scalar kernel regression and `riem.m2skregCV()` for parameter selection using cross validation.

# Riemann 0.1.1

* Interfaces to some of the functions are simplified to emphasize key parameters only for users.
* Added functionalities added *galore*. 
* Vignette for an elementary usage of the package is added.

# Riemann 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Initial release : **Hello, World!**
