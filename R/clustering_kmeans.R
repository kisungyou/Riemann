#' K-Means Clustering
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform k-means clustering by minimizing within-cluster sum of squares (WCSS). 
#' Since the problem is NP-hard and sensitive to the initialization, we provide an 
#' option with multiple starts and return the best result with respect to WCSS.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{algorithm}{(case-insensitive) name of an algorithm; \code{"Lloyd"} (default), or \code{"MacQueen"}.}
#' \item{init}{\code{"plus"} for squared-distance k-means++ (default), \code{"random"} for sampled observations, or a vector of k distinct observation indices.}
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{nstart}{the number of starts (default: 5).}
#' \item{eps}{relative objective-increase tolerance (default: 1e-8).}
#' \item{mean.maxiter}{maximum iterations for each cluster mean (default: 200).}
#' \item{mean.eps}{cluster-mean convergence tolerance (default: 1e-8).}
#' }
#' 
#' @details The default is now Lloyd iteration, with each cluster updated by
#'   the shared Fr\ifelse{html}{\out{&eacute;}}{\ifelse{latex}{\out{\'e}}{e}}chet-mean solver. A Lloyd objective increase beyond tolerance
#'   is rejected and reported. The retained MacQueen option sequentially assigns
#'   observations and recomputes affected cluster means; no global monotonicity
#'   guarantee is claimed for approximate or nonconvex manifold means.
#'   Tied distances choose the lowest cluster index. Empty clusters are reseeded
#'   using a maximal-residual observation from a nonsingleton cluster; impossible
#'   separation and failed mean solves invalidate that start. Returned labels
#'   always assign observations to the final stored centers. Nonconverged finite
#'   starts remain visible, and returning one produces a warning. Use
#'   \code{set.seed()} to reproduce random initializations.
#'
#' @return a fitted \code{riem_kmeans} object containing\describe{
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
#' \item{means}{a 3d array where each slice along 3rd dimension is a matrix representation of class mean.}
#' \item{score}{within-cluster sum of squares (WCSS).}
#' }
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
#' ## K-MEANS WITH K=2,3,4
#' clust2 = riem.kmeans(myriem, k=2)
#' clust3 = riem.kmeans(myriem, k=3)
#' clust4 = riem.kmeans(myriem, k=4)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(mds2d, pch=19, main="true label", col=mylabs)
#' plot(mds2d, pch=19, main="K=2", col=clust2$cluster)
#' plot(mds2d, pch=19, main="K=3", col=clust3$cluster)
#' plot(mds2d, pch=19, main="K=4", col=clust4$cluster)
#' par(opar)
#' 
#' @seealso \code{\link{riem.kmeanspp}}
#' 
#' @references 
#' \insertRef{lloyd_least_1982}{Riemann}
#' 
#' \insertRef{macqueen_methods_1967}{Riemann}
#' 
#' @concept clustering
#' @export
riem.kmeans <- function(riemobj, k = 2, geometry = NULL, ...) {
  riem_validate_data(riemobj)
  geom <- riem_resolve_geometry(riemobj, geometry, capability = "mean")
  n <- length(riemobj$data)
  k <- riem_kmeans_integer(k, "k", maximum = n)
  pars <- list(...)
  allowed <- c("algorithm", "init", "maxiter", "nstart", "eps", "mean.maxiter", "mean.eps")
  if (length(pars) && (is.null(names(pars)) || any(names(pars) == "") ||
                      anyDuplicated(names(pars)) || any(!names(pars) %in% allowed))) {
    stop("Unknown, unnamed, or duplicate clustering control.", call. = FALSE)
  }
  control <- utils::modifyList(list(algorithm = "lloyd", init = "plus", maxiter = 50L,
                                   nstart = 5L, eps = 1e-8, mean.maxiter = 200L,
                                   mean.eps = 1e-8), pars)
  control$algorithm <- match.arg(tolower(control$algorithm), c("lloyd", "macqueen"))
  if (is.numeric(control$init)) {
    if (length(control$init) != k || any(!is.finite(control$init)) ||
        any(control$init != floor(control$init)) || any(control$init < 1 | control$init > n) ||
        anyDuplicated(control$init)) {
      stop("Numeric init must contain k distinct valid observation indices.", call. = FALSE)
    }
  } else control$init <- match.arg(tolower(control$init), c("plus", "random"))
  control$maxiter <- riem_kmeans_integer(control$maxiter, "maxiter")
  control$nstart <- riem_kmeans_integer(control$nstart, "nstart")
  control$mean.maxiter <- riem_kmeans_integer(control$mean.maxiter, "mean.maxiter")
  for (field in c("eps", "mean.eps")) {
    value <- control[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(field, " must be positive and finite.", call. = FALSE)
    }
  }
  records <- vector("list", control$nstart)
  states <- vector("list", control$nstart)
  for (start in seq_len(control$nstart)) {
    result <- tryCatch(riem_kmeans_start(riemobj, k, geom, control$init,
      control$algorithm, control$maxiter, control$eps, control$mean.maxiter,
      control$mean.eps), error = function(error) error)
    failed <- inherits(result, "error")
    records[[start]] <- data.frame(start = start, valid = !failed,
      converged = if (failed) FALSE else result$converged,
      termination = if (failed) "failed" else result$termination,
      objective = if (failed) NA_real_ else result$score,
      iterations = if (failed) NA_integer_ else result$iterations,
      empty_clusters = if (failed) NA_integer_ else length(result$empty_clusters),
      error = if (failed) conditionMessage(result) else "", stringsAsFactors = FALSE)
    states[[start]] <- if (failed) list(error = conditionMessage(result)) else result
  }
  diagnostics <- do.call(rbind, records)
  valid <- which(diagnostics$valid & is.finite(diagnostics$objective))
  if (!length(valid)) {
    condition <- structure(list(message = paste0("All clustering starts failed: ",
      paste(unique(diagnostics$error), collapse = "; ")), call = NULL,
      starts = diagnostics), class = c("riem_kmeans_error", "error", "condition"))
    stop(condition)
  }
  best <- valid[which.min(diagnostics$objective[valid])]
  result <- states[[best]]
  if (!result$converged) {
    warning("The selected clustering start did not converge: ", result$termination,
            ". Inspect $starts and $diagnostics.", call. = FALSE)
  }
  means <- array(unlist(result$centers, use.names = FALSE), dim = c(riemobj$size, k))
  output <- list(call = match.call(), cluster = result$cluster, means = means,
    score = result$score, centers = result$centers, geometry = geom, k = k,
    input_template = riem_input_contract(riemobj), controls = control,
    converged = result$converged, termination = result$termination,
    iterations = result$iterations, starts = diagnostics, start_details = states,
    selected.start = best, objective_history = result$objective_history,
    diagnostics = list(empty_clusters = result$empty_clusters,
      empty_repairs = result$empty_repairs, mean_solver = "shared Frechet mean",
      objective = "unweighted within-cluster sum of squared selected-geometry distances"),
    schema_version = 1L, package_version = as.character(utils::packageVersion("Riemann")))
  structure(output, class = c("riem_kmeans", "riemfit"))
}

#' @rdname riem.kmeans
#' @param object a fitted clustering object.
#' @param newdata compatible wrapped observations.
#' @param type return cluster labels (\code{"class"}) or distances to all trained
#'   centers (\code{"distance"}).
#' @method predict riem_kmeans
#' @export
predict.riem_kmeans <- function(object, newdata, type = c("class", "distance"), ...) {
  riem_kmeans_check_fit(object)
  riem_check_newdata(object$input_template, newdata)
  type <- match.arg(type)
  distances <- riem_kmeans_distances(newdata$data, object$centers, object$geometry)
  colnames(distances) <- paste0("Cluster", seq_len(object$k))
  if (type == "distance") distances else riem_kmeans_assign(distances)
}

#' @rdname riem.kmeans
#' @param x a fitted clustering object.
#' @method print riem_kmeans
#' @export
print.riem_kmeans <- function(x, ...) {
  cat("Manifold k-means:", x$k, "clusters on", x$geometry$geometry_id, "\n")
  cat("Within-cluster sum of squares:", format(x$score),
      "; termination:", x$termination, "\n")
  invisible(x)
}

#' @rdname riem.kmeans
#' @method summary riem_kmeans
#' @export
summary.riem_kmeans <- function(object, ...) {
  riem_kmeans_check_fit(object)
  list(geometry = object$geometry, sizes = tabulate(object$cluster, nbins = object$k),
       score = object$score, converged = object$converged,
       termination = object$termination, starts = object$starts,
       selected.start = object$selected.start, controls = object$controls)
}
