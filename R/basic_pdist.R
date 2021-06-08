#' Compute Pairwise Distances for Data
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, compute 
#' pairwise distances.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) in geometry
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
riem.pdist <- function(riemobj, geometry=c("intrinsic","extrinsic"), as.dist=FALSE){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.pdist : input ",DNAME," should be an object of 'riemdata' class."))
  }
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  mydist     = as.logical(as.dist)

  ## COMPUTE
  distmat = basic_pdist(riemobj$name, riemobj$data, mygeometry)
  if (mydist){
    return(stats::as.dist(distmat))
  } else {
    return(distmat)
  }
}
