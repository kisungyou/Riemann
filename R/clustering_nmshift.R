#' Nonlinear Mean Shift
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform clustering of the data based on the nonlinear mean shift algorithm. 
#' Gaussian kernel is used with the bandwidth \eqn{h} as of 
#' \deqn{G(x_i, x_j) \propto \exp \left( - \frac{\rho^2 (x_i,x_j)}{h^2} \right)}
#' where \eqn{\rho(x,y)} is geodesic distance between two points \eqn{x,y\in\mathcal{M}}. 
#' Numerically, some of the limiting points that collapse into the same cluster are 
#' not exact. For such purpose, we require \code{maxk} parameter to search the 
#' optimal number of clusters based on \eqn{k}-medoids clustering algorithm 
#' in conjunction with silhouette criterion.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param h Finite positive bandwidth for the Gaussian weights.
#' @param maxk Integer upper limit for the silhouette search, at least two.
#'   Candidate counts run from two through the smaller of this limit and one
#'   less than the sample size. At least three observations are required.
#' @param maxiter Positive integer maximum number of updates per observation.
#' @param eps Finite positive tolerance for the change in the matrix representation.
#' 
#' @return a named list containing\describe{
#' \item{distance}{an \eqn{(N\times N)} distance between modes corresponding to each data point.}
#' \item{cluster}{a length-\eqn{N} vector of class labels.}
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
#' set.seed(496)
#' ndata  = 10
#' mydata = list()
#' for (i in 1:ndata){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in (ndata+1):(2*ndata)){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in ((2*ndata)+1):(3*ndata)){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' mylabs = rep(c(1,2,3), each=ndata)
#' 
#' ## RUN NONLINEAR MEANSHIFT FOR DIFFERENT 'h' VALUES
#' run1 = riem.nmshift(myriem, maxk=10, h=0.1)
#' run2 = riem.nmshift(myriem, maxk=10, h=1)
#' run3 = riem.nmshift(myriem, maxk=10, h=10)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' plot(mds2d, pch=19, main="label : h=0.1", col=run1$cluster)
#' plot(mds2d, pch=19, main="label : h=1",   col=run2$cluster)
#' plot(mds2d, pch=19, main="label : h=10",  col=run3$cluster)
#' image(run1$distance[,30:1], axes=FALSE, main="distance : h=0.1")
#' image(run2$distance[,30:1], axes=FALSE, main="distance : h=1")
#' image(run3$distance[,30:1], axes=FALSE, main="distance : h=10")
#' par(opar)
#' 
#' @references 
#' \insertRef{subbarao_nonlinear_2009}{Riemann}
#' 
#' @concept clustering
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
riem.nmshift <- function(riemobj, h = 1, maxk = 5, maxiter = 50, eps = 1e-5) {
  riem_resolve_geometry(riemobj, "intrinsic", capability = "mean")
  n <- length(riemobj$data)
  if (n < 3L) stop("Mean-shift grouping requires at least three observations.", call. = FALSE)
  h <- riem_legacy_positive(h, "h")
  maxiter <- riem_regression_integer(maxiter, "maxiter", 1L)
  eps <- riem_legacy_positive(eps, "eps")
  maxk <- min(riem_regression_integer(maxk, "maxk", 2L), n - 1L)

  shifted <- clustering_nmshift(riemobj$name, riemobj$data, h, maxiter, eps)
  select_medoids <- utils::getFromNamespace("hidden_kmedoids_best", "maotai")
  selection <- select_medoids(stats::as.dist(shifted$distance), 2L, maxk)
  labels <- selection$label[, selection$opt.k - 1L]
  list(distance = shifted$distance, cluster = labels)
}
