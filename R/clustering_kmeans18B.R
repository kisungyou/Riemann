#' K-Means Clustering with Lightweight Coreset
#' 
#' Apply weighted Lloyd iterations to an independently sampled lightweight
#' coreset of manifold-valued observations. The Euclidean coreset approximation
#' theorem is not asserted for arbitrary manifold geometries.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param M integer number of independent draws, at least \eqn{k} (default: the larger of \eqn{k} and \eqn{\lceil N/2 \rceil}). Values greater than \eqn{N} are permitted.
#' @param geometry (case-insensitive) name or saved specification of a geometry supporting means. Legacy aliases \code{"intrinsic"} and \code{"extrinsic"} are accepted; see \code{\link{riem.geometry}}.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of weighted Lloyd passes per start (default:50).}
#' \item{nstart}{the number of random starts (default: 5).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
#' \item{means}{a 3d array where each slice along 3rd dimension is a matrix representation of class mean.}
#' \item{score}{unweighted full-data within-cluster sum of squares (WCSS), used to select the best start.}
#' \item{coreset}{the selected start's sampled indices and importance weights.}
#' \item{iterations,converged,termination}{iteration count, convergence flag, and stopping reason for the selected start.}
#' \item{objective_history}{weighted coreset costs, including the initialization.}
#' \item{empty_clusters}{cluster indices unoccupied by the full data.}
#' \item{starts}{per-start validity, full-data score, convergence, and errors.}
#' }
#' 
#' @details Sampling is with replacement, using the same probabilities and
#' weights as \code{riem.coreset18B}. Repeated indices remain separate weighted
#' observations. Each start uses weighted squared-distance initialization and
#' weighted cluster means. Mean solves use at most 200 iterations and tolerance
#' \code{1e-8}; failed mean solves invalidate that start. A Lloyd step that
#' materially increases weighted coreset cost is rejected. No global minimum
#' is guaranteed. Final labels assign every original observation to its nearest
#' stored center, with distance ties choosing the lowest cluster index.
#' A coreset can contain fewer than \eqn{k} distinct locations; empty centers are
#' retained rather than redrawing and conditioning the sampling distribution.
#' A warning identifies an empty cluster or a nonconverged selected start.
#' Use \code{set.seed()} for reproducibility.
#'
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # class 1 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # class 3 : 10 perturbed data points near (0,0,1) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata = list()
#' for (i in 1:10){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 21:30){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' mylabs = rep(c(1,2,3), each=10)
#' 
#' ## TRY DIFFERENT SIZES OF CORESET WITH K=3 FIXED
#' core1 = riem.kmeans18B(myriem, k=3, M=5)
#' core2 = riem.kmeans18B(myriem, k=3, M=10)
#' core3 = riem.kmeans18B(myriem, k=3, M=15)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(mds2d, pch=19, main="true label", col=mylabs)
#' plot(mds2d, pch=19, main="kmeans18B: M=5",  col=core1$cluster)
#' plot(mds2d, pch=19, main="kmeans18B: M=10", col=core2$cluster)
#' plot(mds2d, pch=19, main="kmeans18B: M=15", col=core3$cluster)
#' par(opar)
#' 
#' @references 
#' \insertRef{bachem_scalable_2018a}{Riemann}
#' 
#' @seealso \code{\link{riem.coreset18B}}
#' @concept clustering
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
riem.kmeans18B <- function(riemobj, k=2,
                          M=max(k, ceiling(length(riemobj$data)/2)),
                          geometry=c("intrinsic","extrinsic"), ...) {
  riem_validate_data(riemobj)
  if (missing(geometry)) geometry <- NULL
  geometry <- riem_resolve_geometry(riemobj, geometry, capability = "mean")
  k <- riem_kmeans_integer(k, "k", maximum = length(riemobj$data))
  M <- riem_kmeans_integer(M, "M", minimum = k)
  pars <- riem_legacy_parameters(list(...), c("maxiter", "nstart"))
  maxiter <- if (is.null(pars$maxiter)) 50L else riem_kmeans_integer(pars$maxiter, "maxiter")
  nstart <- if (is.null(pars$nstart)) 5L else riem_kmeans_integer(pars$nstart, "nstart")
  results <- vector("list", nstart)
  records <- vector("list", nstart)
  for (i in seq_len(nstart)) {
    result <- tryCatch(clustering_kmeans18B(riemobj$name, geometry$backend,
      riemobj$data, k, M, maxiter), error = function(error) error)
    valid <- !inherits(result, "error")
    results[[i]] <- result
    records[[i]] <- data.frame(start = i, valid = valid,
      score = if (valid) result$wcss else NA_real_,
      converged = if (valid) result$converged else FALSE,
      error = if (valid) "" else conditionMessage(result), stringsAsFactors = FALSE)
  }
  starts <- do.call(rbind, records)
  valid <- which(starts$valid & is.finite(starts$score))
  if (!length(valid)) stop("All coreset clustering starts failed: ",
    paste(unique(starts$error), collapse = "; "), call. = FALSE)
  best <- results[[valid[which.min(starts$score[valid])]]]
  cluster <- as.integer(best$cluster) + 1L
  empty <- which(tabulate(cluster, nbins = k) == 0L)
  if (!isTRUE(best$converged)) warning("The selected coreset clustering start did not converge: ",
    best$termination, ".", call. = FALSE)
  if (length(empty)) warning("The selected coreset leaves ", length(empty),
    " empty cluster(s); increase M or reduce k if distinct clusters are required.", call. = FALSE)
  list(cluster = cluster, means = best$means, score = as.double(best$wcss),
    coreset = list(coreid = as.integer(best$coreid) + 1L, weight = as.vector(best$weight)),
    iterations = best$iterations, converged = best$converged, termination = best$termination,
    objective_history = as.vector(best$objective_history), empty_clusters = empty,
    starts = starts)
}
