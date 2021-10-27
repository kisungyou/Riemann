#' Pairwise Distance on SPD Manifold
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N} in SPD manifold, compute 
#' pairwise distances among observations.
#'
#' @param spdobj a S3 \code{"riemdata"} class of SPD-valued data.
#' @param geometry name of the geometry to be used. See \code{\link{spd.geometry}} for supported geometries.
#' @param as.dist logical; if \code{TRUE}, it returns a \code{dist} object. Else, it returns a symmetric matrix.
#'
#' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #                   Two Types of Covariances
#' #
#' #  group1 : perturbed from data by N(0,1) in R^3
#' #  group2 : perturbed from data by [sin(x); cos(x); sin(x)*cos(x)]
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' spd_mats = array(0,c(3,3,20))
#' for (i in 1:10){
#'   spd_mats[,,i] = stats::cov(matrix(rnorm(50*3), ncol=3))
#' }
#' for (j in 11:20){
#'   randvec = stats::rnorm(50, sd=3)
#'   randmat = cbind(sin(randvec), cos(randvec), sin(randvec)*cos(randvec))
#'   spd_mats[,,j] = stats::cov(randmat + matrix(rnorm(50*3, sd=0.1), ncol=3))
#' }
#' 
#' ## WRAP IT AS SPD OBJECT
#' spd_obj = wrap.spd(spd_mats)
#' 
#' ## COMPUTE PAIRWISE DISTANCES
#' #  Geometries are case-insensitive.
#' pdA = spd.pdist(spd_obj, "airM")
#' pdL = spd.pdist(spd_obj, "lErm")
#' pdJ = spd.pdist(spd_obj, "Jeffrey")
#' pdS = spd.pdist(spd_obj, "stEin")
#' pdW = spd.pdist(spd_obj, "wasserstein")
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' image(pdA, axes=FALSE, main="AIRM")
#' image(pdL, axes=FALSE, main="LERM")
#' image(pdJ, axes=FALSE, main="Jeffrey")
#' image(pdS, axes=FALSE, main="Stein")
#' image(pdW, axes=FALSE, main="Wasserstein")
#' par(opar)
#' }
#'
#' @concept spd
#' @export
spd.pdist <- function(spdobj, geometry, as.dist=FALSE){
    # PREPARE
    DNAME = paste0("'",deparse(substitute(spdobj)),"'")
    if ((!inherits(spdobj,"riemdata"))||(!all(spdobj$name=="spd"))){
      stop(paste0("* spd.pdist : input ",DNAME," should be an object of 'riemdata' class on 'spd' manifold.."))
    }
    mygeom = spd.geometry(geometry)
    mydist = as.logical(as.dist)
    
    # COMPUTE
    array3d = spd.wrap3d(spdobj$data)
    output  = src_spd_pdist(array3d, mygeom)
    
    # RETURN
    if (mydist){
      return(stats::as.dist(output))
    } else {
      return(output)
    }
}