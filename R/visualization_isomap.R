#' Isometric Feature Mapping
#' 
#' ISOMAP - isometric feature mapping - is a dimensionality reduction method 
#' to apply classical multidimensional scaling to the geodesic distance 
#' that is computed on a weighted nearest neighborhood graph. Nearest neighbor 
#' is defined by \eqn{k}-NN where two observations are said to be connected when 
#' they are mutually included in each other's nearest neighbor. Note that 
#' it is possible for geodesic distances to be \code{Inf} when nearest neighbor 
#' graph construction incurs separate connected components. When an extra 
#' parameter \code{padding=TRUE}, infinite distances are explicitly replaced by 2 times
#' the maximal finite geodesic distance.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param nnbd the size of nearest neighborhood (default: 5).
#' @param geometry A geometry name or saved specification.
#' @details Self-neighbors are excluded and distance ties use observation order.
#'   An edge is retained only when each endpoint selects the other among its
#'   nearest neighbors. Zero-length edges between duplicate observations are
#'   retained. The result records components, adjacency, and whether padding
#'   changed disconnected distances. Padding produces an artificial dissimilarity
#'   and is not an estimate of the original disconnected graph distance.
#' @param ... extra parameters including\describe{
#' \item{padding}{a logical, default \code{FALSE}; if \code{TRUE}, disconnected
#' distances are replaced by twice the largest finite graph distance, with a warning.}
#' }
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # 10 perturbed data points near (0,0,1) on S^2 in R^3
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
#' ## MDS AND ISOMAP WITH DIFFERENT NEIGHBORHOOD SIZE
#' mdss = riem.mds(myriem)$embed
#' iso1 = riem.isomap(myriem, nnbd=5, padding=TRUE)$embed
#' iso2 = riem.isomap(myriem, nnbd=10, padding=TRUE)$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(mdss, col=mylabs, pch=19, main="MDS")
#' plot(iso1, col=mylabs, pch=19, main="ISOMAP:nnbd=5")
#' plot(iso2, col=mylabs, pch=19, main="ISOMAP:nnbd=10")
#' par(opar)
#' 
#' @references
#' Silva VD and Tenenbaum JB (2003). "Global Versus Local Methods in Nonlinear
#' Dimensionality Reduction." Advances in Neural Information Processing Systems
#' 15, 721--728. MIT Press.
#' 
#' @concept visualization
#' @export
riem.isomap <- function(riemobj, ndim = 2, nnbd = 5, geometry = NULL, ...) {
  result <- riem_legacy_distances(riemobj, geometry)
  n <- length(riemobj$data)
  ndim <- riem_legacy_dimensions(ndim, n)
  nnbd <- riem_regression_integer(nnbd, "nnbd", 1L, n - 1L)
  options <- list(...)
  if (length(options) && (is.null(names(options)) || any(!names(options) %in% "padding") ||
                          anyDuplicated(names(options)))) {
    stop("Only the named 'padding' option is supported.", call. = FALSE)
  }
  padding <- if (is.null(options$padding)) FALSE else options$padding
  if (!is.logical(padding) || length(padding) != 1L || is.na(padding)) {
    stop("padding must be TRUE or FALSE.", call. = FALSE)
  }
  graph <- riem_mutual_knn_graph(result$distances, nnbd)
  disconnected <- any(!is.finite(graph$distance))
  if (disconnected) {
    if (!padding) stop("The mutual-neighbor graph is disconnected; increase nnbd or explicitly request padding=TRUE.",
                       call. = FALSE)
    fill <- 2 * max(graph$distance[is.finite(graph$distance)])
    if (!is.finite(fill)) stop("The requested padding distance is not finite.", call. = FALSE)
    graph$distance[!is.finite(graph$distance)] <- fill
    warning("Disconnected graph distances were padded; the embedding uses an artificial dissimilarity.", call. = FALSE)
  }
  out <- riem_legacy_cmds(graph$distance, ndim)
  out$geometry <- result$geometry
  out$adjacency <- graph$adjacency
  out$component <- graph$component
  out$padded <- disconnected
  out$nnbd <- nnbd
  out$method <- "mutual_knn_isomap"
  out
}
