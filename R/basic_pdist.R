#' Compute Pairwise Distances for Data
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, compute 
#' pairwise distances.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry A name or saved specification accepted by \code{\link{riem.geometry}}.
#'   NULL uses the object's saved geometry, or intrinsic geometry for an ordinary wrapper.
#' @param as.dist logical; if \code{TRUE}, it returns \code{dist} object, else it returns a symmetric matrix.
#' 
#' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with two types
#' #
#' #  group1 : perturbed data points near (0,0,1) on S^2 in R^3
#' #  group2 : perturbed data points near (1,0,0) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata = list()
#' sdval  = 0.1
#' for (i in 1:10){
#'   tgt = c(stats::rnorm(2, sd=sdval), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(1, stats::rnorm(2, sd=sdval))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' 
#' ## COMPARE TWO DISTANCES
#' dint = riem.pdist(myriem, geometry="intrinsic", as.dist=FALSE)
#' dext = riem.pdist(myriem, geometry="extrinsic", as.dist=FALSE)
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dint[,nrow(dint):1], main="intrinsic", axes=FALSE)
#' image(dext[,nrow(dext):1], main="extrinsic", axes=FALSE)
#' par(opar)
#' 
#' @concept basic
#' @export
riem.pdist <- function(riemobj, geometry = NULL, as.dist = FALSE) {
  spec <- riem_resolve_geometry(riemobj, geometry, "distance")
  if (!is.logical(as.dist) || length(as.dist) != 1L || is.na(as.dist)) {
    stop("as.dist must be TRUE or FALSE.", call. = FALSE)
  }
  out <- basic_pdist(riemobj$name, riemobj$data, spec$backend)
  if (any(!is.finite(out)) || any(out < 0)) stop("Distance computation returned invalid values.", call. = FALSE)
  if (as.dist) stats::as.dist(out) else out
}
