#!/usr/bin/env Rscript
# Explicit per-export disposition. No package functions are evaluated.
# An inventory run supplies its resolved context explicitly through sys.source().
if (exists('riemann_contract_context', inherits = FALSE)) {
  maintenance_context <- riemann_contract_context
} else {
  maintenance_script <- local({
    sourced <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
    if (length(sourced)) return(tail(sourced, 1L))
    invoked <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(invoked)) stop('Run this entry point with Rscript or source().')
    sub('^--file=', '', invoked[1L])
  })
  source(file.path(dirname(normalizePath(maintenance_script, mustWork = TRUE)), 'paths.R'),
         local = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  maintenance_context <- riemann_maintenance_paths(maintenance_script,
    if (length(args)) args[1L] else NULL, require_package = TRUE)
}
contract_root <- maintenance_context$package_root
contract_inventory <- maintenance_context$inventory
if (is.null(contract_inventory)) {
  stop('Run tools/maintenance/generate-inventory.R to regenerate method contracts.')
}
contract_rows <- list()
record_contract <- function(names, disposition, domain, assumptions, limitations, evidence) {
  for (name in names) contract_rows[[name]] <<- data.frame(
    name = name, disposition = disposition, domain = domain,
    mathematical_or_statistical_assumptions = assumptions,
    limitations = limitations, evidence_scope = evidence, stringsAsFactors = FALSE)
}
record_contract(grep('^wrap\\.', contract_inventory$name, value = TRUE), 'core',
  'Input representation only; accepted dimensions, finite values and manifold membership are validated.',
  'A valid wrapped representation does not establish availability or correctness of a downstream geometry or estimator.',
  'Correlation distance/summary, intrinsic Stiefel distance/summary, extrinsic fixed-rank SPD and extrinsic landmark summaries are restricted separately by the geometry layer.',
  'Wrapper validation source and test-wrappers/test-membership-persistence evidence; consult current package test logs.')
record_contract(c('riem.geometry', 'riem.capabilities', 'print.riem_geometry'), 'core',
  'Serializable geometry metadata and primitive availability registry.',
  'Metric, representation dimension and persisted parameters belong to the statistical model.',
  'Availability is distinct from mathematical validation and method calibration; registry validation_status must be read.',
  'Geometry specification and persistence tests.')
record_contract(c('riem.pdist', 'riem.pdist2'), 'core',
  'Available distance backends resolved for compatible riemdata.',
  'Distances are unsquared; quotient identification and invariances follow the resolved geometry.',
  'Core independent geometry checks emphasize Euclidean, SPD and sphere; other available native manifolds retain their recorded validation status.',
  'Primitive boundary/independent-reference tests; full manifold certification is not inferred.')
record_contract(c('riem.mean', 'riem.median'), 'core',
  'Resolved mean/median capabilities; finite nonnegative weights of positive total mass.',
  'Objective uses the same declared distances. Intrinsic solvers return numerical estimates and convergence diagnostics; stationary termination establishes a first-order condition only.',
  'A stationary point can be a saddle or maximum; neither local minimality nor uniqueness/global optimality is certified. Some extrinsic medians are explicitly projected ambient estimates. Unsupported geometry combinations error.',
  'Shared solver tests and low-dimensional/reference fixtures.')
record_contract(c('print.riemdata', 'summary.riemdata', 'print.riem_summary',
                   'summary.riem_summary', 'plot.riem_summary'), 'core',
  'Presentation of validated inputs and saved summary metadata.',
  'Presentation adds no inferential calibration.',
  'Convergence and estimand fields of the underlying fit remain authoritative.',
  'Summary object/interface tests; no separate estimator claim.')
record_contract(grep('m2skreg', contract_inventory$name, value = TRUE), 'core',
  'Scalar-response Gaussian Nadaraya-Watson regression using available distances.',
  'Finite scalar responses; positive finite bandwidth; explicit fold assignment or randomized K-fold validation.',
  'Self-weights enter training fitted values. Positive radial weights do not assert that a kernel is positive semidefinite. Prediction stores geometry; no automatic inferential intervals.',
  'Independent weighted-average and all-fold CV fixtures, extreme bandwidth/response tests, failure and persistence tests.')
record_contract(c('riem.kmeans', 'riem.kmeanspp', 'predict.riem_kmeans',
                   'summary.riem_kmeans', 'print.riem_kmeans'), 'core',
  'Squared-distance clustering with shared compatible means and explicit algorithms.',
  'Lloyd/local optimization, initialization, ties and empty clusters follow the saved controls.',
  'General-manifold global optimum and approximation guarantees are not claimed; accelerated pruning requires the documented metric contract.',
  'Clustering objective, prediction and diagnostic tests; read fit convergence.')
record_contract(c('riem.pga', 'predict.riem_pga', 'riem.reconstruct',
                   'riem.reconstruct.riem_pga', 'print.riem_pga', 'summary.riem_pga',
                   'print.summary.riem_pga', 'plot.riem_pga'), 'core',
  'Declared tangent/chart PCA: Euclidean, SPD affine-invariant or log-Euclidean, round sphere and supported regular landmark configurations.',
  'Coordinates are metric compatible at the saved base point; prediction reuses that base point and basis.',
  'Tangent/chart approximation rather than exact nonlinear PGA. Cut loci and singular/ambiguous landmark alignments are restricted; explained variance is in the declared coordinate model.',
  'Independent coordinate, reconstruction, variance and persistence tests, including highly unequal positive weights and normalization underflow rejection.')
record_contract('riem.mds', 'validated_restricted', 'Classical scaling of a finite available distance matrix.',
  'Double-centering squared distances; descending positive eigenspace with sqrt-eigenvalue scores.',
  'Non-Euclidean distances can have a negative spectrum, which is returned and explicitly truncated or rejected; no exact isometry claim. Requested excess dimensions are zero-filled.',
  'Euclidean reconstruction, constant/rank-deficient input and spherical negative-spectrum tests.')
record_contract('riem.kpca', 'validated_restricted', 'Centered Gaussian similarity based on the legacy extrinsic distance definition.',
  'A reproducing-kernel interpretation requires a positive-semidefinite kernel, not merely nonnegative pair weights.',
  'Empirical negative spectrum is reported, warns and is explicitly truncated or rejected. Truncation is an approximation; an empirically PSD sample does not prove global positive definiteness.',
  'Independent Gram reconstruction, sqrt-eigenvalue scores, spectrum ordering and indefinite-matrix tests.')
record_contract('riem.isomap', 'validated_restricted', 'Classical scaling of shortest paths on an actual mutual k-nearest-neighbor graph.',
  'Graph approximation depends on sample coverage and neighborhood choice; zero-cost duplicate edges are retained.',
  'Disconnected graphs error unless artificial padding is explicitly requested and reported. No universal geodesic recovery or global isometry guarantee.',
  'Independent small-graph adjacency/components, duplicate-edge and complete-graph tests.')
record_contract('riem.hclust', 'validated_restricted', 'stats::hclust applied to available unsquared dissimilarities.',
  'Linkage-specific interpretation follows stats::hclust.',
  'Ward, centroid and median Euclidean variance interpretations do not transfer automatically to arbitrary manifold distances; no Frechet-centroid claim.',
  'Supported stats interface comparison; validated finite distances and members.')
record_contract(c('riem.fanova', 'riem.fanovaP'), 'validated_restricted',
  'Compatible converged Frechet means and nondegenerate empirical squared-distance variances, at least three observations per group.',
  'Asymptotic test compares Frechet means and variances under uniqueness, moments/entropy and nonvanishing group proportions. Label permutation requires independent exchangeable observations under identical group distributions.',
  'Not an omnibus all-distribution test; paired/clustered data require another resampling design. Empirical nondegeneracy and solver convergence do not prove population assumptions.',
  'Independent Euclidean and commuting SPD log-coordinate statistics; proportions, scale invariance, refit permutations and plus-one p-value fixtures.')
record_contract('riem.test2bg14', 'validated_restricted',
  'Permutation statistic of within/between mean distances, at least two observations per group.',
  'Independent observations exchangeable under identical group distributions; fixed geometry for every assignment.',
  'The tested mean-distance discrepancy is not asserted to detect every alternative in every manifold. No finite-population or dependent-data calibration.',
  'Independent statistic, compatible groups, inclusive ties and plus-one permutation tests.')
record_contract(c('riem.wasserstein', 'riem.test2wass'), 'validated_restricted',
  'Empirical probability measures with finite available ground distances and finite p >= 1.',
  'Exact default calls the external linear-program solver and checks plan feasibility; test permutations require independent exchangeable observations.',
  'Numerical plan checks do not certify optimality. Explicit use.smooth=TRUE remains an experimental IPOT approximation with no optimality certificate; high powers can lose precision.',
  'One-dimensional transport reference, extreme weight normalization, singleton matrices, marginals and permutation fixtures.')
record_contract(c('dacg', 'dmacg', 'racg', 'rmacg', 'mle.acg', 'mle.macg'), 'validated_restricted',
  'Angular models on sphere/Stiefel (MACG also descends to Grassmann), densities relative to uniform probability measure.',
  'Positive-definite shape parameter is scale invariant and trace normalized. Interior identifiable estimation requires adequate dispersion across proper subspaces.',
  'Fixed-point convergence is reported but does not certify existence/uniqueness. Full-dimensional MACG estimation is unidentifiable and errors. Density overflow may require log=TRUE.',
  'Independent density formula/integration, scale and basis invariance, polar sample membership and fixed-point residual tests.')
record_contract(c('riem.scUL', 'riem.scSM', 'riem.scNJW', 'riem.sc05Z'), 'experimental',
  'Spectral clustering from nonnegative affinities via T4cluster.',
  'Graph normalization and connectedness/positive degrees matter. Affinity graphs do not require a PSD Mercer kernel.',
  'External implementation, zero-degree/disconnected cases and full parameter boundaries are not comprehensively validated. Partition is descriptive, not inferential evidence.',
  'Static source and dependency audit; broad independent mathematical validation pending.')
record_contract('riem.kmedoids', 'experimental', 'PAM-style dissimilarity clustering via a private maotai helper.',
  'Medoids are observed objects and dissimilarities determine the objective.',
  'Private dependency interface and algorithm diagnostics require consolidation; a global optimum is not promised.',
  'Euclidean fixture checks assignments against distances to returned medoids; independent optimality assessment remains pending.')
record_contract(c('riem.coreset18B', 'riem.kmeans18B'), 'experimental',
  'Importance-sampled weighted coreset and associated clustering.',
  'Independent draws with replacement use weights 1/(M q_i); clustering applies those weights to compatible Frechet center updates. Cited approximation bounds depend on the source theorem distance/center assumptions.',
  'Theorem conditions are not proved for every native manifold. Duplicate draws and empty clusters are possible; a finite local fit is not a global optimization or approximation certificate.',
  'Independent finite-design expectation, sampling weight, weighted-center/objective, zero-dispersion and iteration-control regression tests; theorem transfer remains unestablished.')
record_contract('riem.clrq', 'experimental', 'Competitive learning by stochastic intrinsic updates.',
  'Intrinsic exp/log/metric consistency and a valid stochastic approximation gain schedule are needed.',
  'Default t^-1/2-type gain does not by itself meet the usual square-summable stochastic approximation condition; no general convergence claim. Private initialization helper remains.',
  'Euclidean fixture checks assignments after private-backend initialization; general convergence remains unverified.')
record_contract('riem.nmshift', 'experimental', 'Intrinsic non-blurring mean-shift modes with downstream medoid grouping.',
  'A radial weighting rule is distinct from a normalized manifold density and must be interpreted with its metric.',
  'Global density-mode consistency and bandwidth effects are not established for every manifold; private model-selection helper remains.',
  'Euclidean reference checks the first Gaussian update, the search bound and a narrow bandwidth; statistical mode calibration remains unverified.')
record_contract('riem.knn', 'experimental', 'Neighbor indices/distances using the selected native distance.',
  'Finite compatible observations; 1 <= k < N. Self is excluded by index; distinct zero-distance points are eligible and distance ties use input index.',
  'Distance availability and numerical limitations inherit the selected geometry; no prediction or inferential calibration is provided.',
  'Independent neighbor ordering, duplicated-point/self-exclusion, count, singleton and geometry-boundary tests.')
record_contract('riem.rmml', 'experimental', 'Label-supervised intrinsic distance transformation.',
  'Finite compatible labels and the metric-learning objective/regularizer must justify any reported metric.',
  'Positive semidefiniteness, label-degeneracy behavior and out-of-sample statistical contract require independent review.',
  'Independent Grassmann projector embedding, basis invariance and regularized metric-learning fixtures; broader learned-distance validity not certified.')
record_contract('riem.seb', 'experimental', 'Iterative smallest-enclosing-ball approximation.',
  'Geodesic convexity, curvature and radius assumptions of the source method determine convergence guarantees.',
  'No general-manifold global optimality/certified enclosure claim; solver and boundary review pending.',
  'Static source/native call mapping.')
record_contract(c('riem.sammon', 'riem.tsne', 'riem.phate'), 'experimental',
  'Descriptive nonlinear ordination from pairwise distances.',
  'Sammon requires a policy for zero dissimilarities; t-SNE uses neighborhood probabilities; PHATE uses a diffusion graph.',
  'Low-dimensional visual separation is not statistical group evidence. Parameter/duplicate/disconnected-graph behavior and private helper dependencies remain to be consolidated; no out-of-sample model is supplied.',
  'Fixed-seed t-SNE and PHATE fixtures check finite embeddings from private backends. Independent validation across inputs remains pending.')
record_contract(c('riem.interp', 'riem.interps'), 'experimental',
  'Selected native interpolation between manifold observations.',
  'Geodesic uniqueness and the selected intrinsic/extrinsic interpolation definition must hold.',
  'Cut-locus, endpoint, multi-segment and out-of-range behavior are not comprehensively audited; a wrapper passing validation does not ensure a unique path.',
  'Static source/native map; geometry-level restrictions apply independently.')
record_contract(c('riem.distlp', 'riem.dtw'), 'experimental',
  'Sampled curve discrepancy: numerical Lp integration or dynamic time warping.',
  'Lp requires a common ordered sampling grid; DTW depends on an explicit step pattern and endpoint convention.',
  'Grid/zero-length/singleton behavior and numerical integration need review. DTW is generally a dissimilarity, not a metric, and should not support metric-tree guarantees.',
  'Static source/native review; independent curve fixtures pending.')
record_contract(c('dsplaplace', 'rsplaplace', 'mle.splaplace', 'dspnorm', 'rspnorm', 'mle.spnorm'),
  'experimental', 'Geodesic spherical Laplace or spherical normal distributions.',
  'Sphere dimension, reference surface measure, normalizer and positive scale/concentration determine the density.',
  'Normalizer accuracy, extreme concentration, sampler distribution and optimizer convergence are not yet comprehensively independently verified.',
  'Static source and native dependency mapping; full distributional validation pending.')
record_contract(c('moSL', 'moSN', 'density', 'density.moSL', 'density.moSN', 'label',
                   'label.moSL', 'label.moSN', 'loglkd', 'loglkd.moSL', 'loglkd.moSN'),
  'experimental', 'Spherical Laplace/normal mixture EM and associated accessors.',
  'Component normalization, finite likelihood, mixing weights and local EM optimization must be validated.',
  'Mixture identifiability, local optima, likelihood monotonicity, concentration boundaries and information-criterion parameter counts need audit. Accessors inherit the fitted-model limitations.',
  'Static source/dependency mapping; no universal mixture fitting or selection validation.')
record_contract(c('sphere.runif', 'stiefel.runif', 'grassmann.runif'), 'experimental',
  'Uniform probability draws in the declared representation.',
  'Sphere normalization and Stiefel/Grassmann orthogonalization must generate the intended invariant probability law.',
  'Representation validity alone does not verify a distribution or cross-platform RNG sequence; full independent moment/distribution tests remain pending.',
  'Source/native RNG mapping and manifold membership checks where available.')
record_contract(c('sphere.utest', 'stiefel.utest', 'grassmann.utest'), 'experimental',
  'Asymptotic Rayleigh or Bingham moment tests of invariant uniformity.',
  'Independent identically distributed observations, appropriate dimensions and asymptotic calibration; moment tests need not detect all nonuniform alternatives.',
  'Finite-sample and modified-expansion calibration are not certified. Grassmann proper-subspace requirement and singular p=2 modified formula are explicitly restricted.',
  'Static formula review plus explicit singular-domain boundary tests; Monte Carlo size/power study pending.')
record_contract(c('sphere.geo2xyz', 'sphere.xyz2geo'), 'experimental',
  'Geographic/cartesian coordinate conversion on the unit two-sphere.',
  'Documented longitude/latitude ordering, decimal degrees and unit-vector inputs.',
  'Pole/antimeridian, invalid norm, recycling and input-shape contracts remain independently unaudited; this is representation conversion rather than inference.',
  'Static source mapping; systematic conversion fixtures pending.')
record_contract(c('grassmann.optmacg', 'stiefel.optSA'), 'experimental',
  'Stochastic black-box optimization on Grassmann/Stiefel representations.',
  'Objective must be finite scalar and, for Grassmann, invariant to the chosen orthonormal basis.',
  'No global optimum certificate. Population covariance fitting/annealing controls and stochastic stopping require further independent review.',
  'Static algorithm/dependency mapping; optimization guarantees pending.')
record_contract(c('spd.geometry', 'spd.pdist', 'spd.wassbary'), 'experimental',
  'Legacy named SPD distances or Bures-Wasserstein barycenter interface.',
  'Each named dissimilarity/metric and barycenter has its own mathematical definition; names do not establish equivalence.',
  'These legacy choices are distinct from the consolidated saved affine-invariant/log-Euclidean models. Metric properties, numerical stability and barycenter convergence require individual audits.',
  'Static source/native mapping; broad legacy metric validation pending.')
record_contract('rmvnorm', 'experimental', 'Gaussian simulation utility.',
  'Finite dimension-compatible mean and positive-definite covariance with positive integer sample count.',
  'Legacy boundary validation and single-versus-multiple return shape need consolidation. Native covariance/RNG correctness is not inferred from compilation.',
  'Static source/native mapping; independent distribution/boundary audit pending.')
# Narrow, dated numerical evidence supplements the broader legacy dispositions.
record_contract('riem.sammon', 'experimental',
  'Descriptive sample-only Sammon ordination from finite pairwise distances.',
  'Off-diagonal original distances must be strictly positive. Optimization uses normalized distances, the Sammon diagonal Hessian and descent backtracking.',
  'No global minimum or class-recovery guarantee. Duplicate observations error because the inverse-distance loss is undefined; no out-of-sample model is supplied.',
  'Independent loss and derivative calculations, scale and descent checks, duplicate policy and the October 2026 divergence regression fixture.')
record_contract(c('dspnorm', 'dsplaplace', 'mle.spnorm', 'mle.splaplace'), 'experimental',
  'Geodesic spherical normal/Laplace likelihoods relative to surface-area measure.',
  'Stable radial quadrature and direct log densities; dispersion estimation solves the score on an adaptive bracket with explicit uniform and point-mass boundaries.',
  'Location fitting on a sphere need not find a global optimum. Point-mass limits are not ordinary surface densities. Broad high-dimensional optimizer/distributional certification is not claimed.',
  'Independent circle normalizers, likelihood scores, concentration/scale boundaries, finite log-tail and optimizer agreement regression tests.')
record_contract(c('grassmann.optmacg', 'stiefel.optSA'), 'experimental',
  'Stochastic black-box optimization on Grassmann/Stiefel representations.',
  'Objective must be a finite scalar and Grassmann objectives must be basis invariant. Stiefel QR factors receive column signs from the triangular diagonal.',
  'No global optimum certificate; finite-run annealing/population convergence is not established.',
  'Tests include both Stiefel hemispheres, manifold membership and the attained cost of a nonconstant Grassmann objective with a known optimum.')
record_contract('spd.pdist', 'experimental',
  'Legacy named SPD distances, with affine-invariant/log-Euclidean choices routed to the shared distance interface.',
  'Stein uses the scaled relative spectrum; Bures-Wasserstein uses a Cholesky-factor orthogonal-alignment residual. Inputs must remain numerically positive definite.',
  'Extreme conditioning can prevent factorization. Named legacy dissimilarities do not define interchangeable fitted-model geometries, and not every choice has received the same validation.',
  'Independent diagonal/2-by-2 formulas, scaling/congruence checks, identical/nearby matrix and determinant-overflow regression tests for Stein and Bures-Wasserstein.')
missing_contracts <- setdiff(contract_inventory$name, names(contract_rows))
if (length(missing_contracts)) stop('Explicit disposition missing: ', paste(missing_contracts, collapse = ', '))
contract_table <- do.call(rbind, contract_rows[contract_inventory$name])
rownames(contract_table) <- NULL
contract_table$source <- contract_inventory$source
contract_table$direct_test_files <- contract_inventory$direct_test_files
contract_table$geometry_status_separate <- TRUE
stopifnot(!anyDuplicated(contract_table$name), all(contract_table$disposition %in%
  c('core', 'validated_restricted', 'experimental', 'disabled')))
dir.create(file.path(contract_root, 'inst'), showWarnings = FALSE)
write.csv(contract_table, file.path(contract_root, 'inst', 'method-contracts.csv'), row.names = FALSE)
print(table(contract_table$disposition))
