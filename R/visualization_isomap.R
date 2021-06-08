#' Isometric Feature Mapping
#' 
#' ISOMAP - isometric feature mapping - is a dimensionality reduction method 
#' to apply classical multidimensional scaling to the geodesic distance 
#' that is computed on a weighted nearest neighborhood graph. Nearest neighbor 
#' is defined by \eqn{k}-NN where two observations are said to be connected when 
#' they are mutually included in each other's nearest neighbor. Note that 
#' it is possible for geodesic distances to be \code{Inf} when nearest neighbor 
#' graph construction incurs separate connected components. When an extra 
#' parameter \code{padding=TRUE}, infinite distances are replaced by 2 times 
#' the maximal finite geodesic distance.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param nnbd the size of nearest neighborhood (default: 5).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{padding}{a logical; if \code{TRUE}, \code{Inf}-valued geodesic distances are replaced by 2 times the maximal geodesic distance in the data.}
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
#' iso1 = riem.isomap(myriem, nnbd=5)$embed
#' iso2 = riem.isomap(myriem, nnbd=10)$embed
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
#' \insertRef{silva_global_2003}{Rdimtools}
#' 
#' @concept visualization
#' @export
riem.isomap <- function(riemobj, ndim=2, nnbd=5, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.mds : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim = max(2, round(ndim))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  mynnbd = max(2, round(nnbd))
  
  ## IMPLICIT PARAMETERS
  params = list(...)
  pnames = names(params)
  use.padding = ifelse(("padding"%in%pnames), as.logical(params$padding), TRUE)
  
  ## COMPUTE WEIGHTED PAIRWISE DISTANCE
  distobj = stats::as.dist(visualize_isomap(riemobj$name, riemobj$data, mygeom, mynnbd))
  # distgeo = maotai::shortestpath(distobj)
  distgeo = Rdimtools::aux.shortestpath(distobj)
  if (any(is.infinite(distgeo))){
    if (use.padding){
      print("* riem.isomap : some of the geodesic distances are Inf, so 'padding' is applied.")  
      distgeo[is.infinite(distgeo)] = max(distgeo[!is.infinite(distgeo)])*2
    } else {
      stop("* riem.isomap : some of the points are isolated. Use larger 'nnbd' value.")
    }
  }
  
  
  ## COMPUTE MDS AND RETURN
  func.import     = utils::getFromNamespace("hidden_cmds", "maotai")
  out.cmds        = func.import(stats::as.dist(distgeo), ndim=myndim)
  out.cmds$stress = NULL
  return(out.cmds)
}
