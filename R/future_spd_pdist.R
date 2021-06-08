#' #' Pairwise Distance on SPD Manifold
#' #' 
#' #' 
#' #' 
#' #' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' #' 
#' #' @concept spd
#' #' @keywords internal
#' #' @noRd
#' spd.pdist <- function(riemobj, geometry, as.dist=FALSE){
#'   # PREPARE
#'   DNAME = paste0("'",deparse(substitute(riemobj)),"'")
#'   if ((!inherits(riemobj,"riemdata"))||(!all(riemobj$name=="spd"))){
#'     stop(paste0("* spd.pdist : input ",DNAME," should be an object of 'riemdata' class on 'spd' manifold.."))
#'   }
#'   mygeom = spd.geometry(geometry)
#'   mydist = as.logical(as.dist)
#'   
#'   # COMPUTE
#'   distmat = spd_pdist(riemobj$data, mygeom)
#'   if (mydist){
#'     return(stats::as.dist(distmat))
#'   } else {
#'     return(distmat)
#'   }
#' }
