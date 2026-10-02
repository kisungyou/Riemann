#' Hierarchical Agglomerative Clustering
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}}, 
#' perform hierarchical agglomerative clustering with 
#' \code{stats::hclust}. The supplied manifold distances are treated as dissimilarities.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry A geometry name or saved specification.
#' @details Ward, centroid and median linkage require care when interpreted as
#'   Euclidean sums of squares. In particular, Ward linkage of geodesic distances
#'   is not manifold k-means. \code{ward.D} is the historical R update and does
#'   not implement the Ward criterion implemented by \code{ward.D2}. See the
#'   \code{stats::hclust} documentation for distance powers and member weights.
#' @param method agglomeration method to be used. This must be one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
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
#' 
#' ## COMPUTE SINGLE AND COMPLETE LINKAGE
#' hc.sing <- riem.hclust(myriem, method="single")
#' hc.comp <- riem.hclust(myriem, method="complete")
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(hc.sing, main="single linkage")
#' plot(hc.comp, main="complete linkage")
#' par(opar)
#' 
#' @references 
#' \insertRef{mullner_fastcluster_2013}{Riemann}
#' 
#' @concept clustering
#' @export
riem.hclust <- function(riemobj, geometry = NULL,
                        method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                   "centroid", "median"), members = NULL) {
  result <- riem_legacy_distances(riemobj, geometry)
  if (length(riemobj$data) < 2L) stop("Hierarchical clustering requires at least two observations.", call. = FALSE)
  method <- match.arg(method)
  if (!is.null(members) && (!is.numeric(members) || is.complex(members) ||
      length(members) != length(riemobj$data) || any(!is.finite(members)) || any(members <= 0))) {
    stop("members must be finite positive cluster sizes, one per input.", call. = FALSE)
  }
  out <- stats::hclust(stats::as.dist(result$distances), method = method, members = members)
  out$geometry <- result$geometry
  out$distance_interpretation <- "supplied dissimilarity; no manifold variance-minimization claim"
  out$call <- match.call()
  out
}
