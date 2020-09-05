#' Multidimensional Scaling
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' apply multidimensional scaling to get low-dimensional embedding 
#' in Euclidean space. Usually, \code{ndim=2,3} are chosen for visualization.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' \item{stress}{discrepancy between embedded and original distances as a measure of error.}
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
#' ## MDS EMBEDDING WITH TWO GEOMETRIES
#' embed2int = riem.mds(myriem, geometry="intrinsic")$embed
#' embed2ext = riem.mds(myriem, geometry="extrinsic")$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2int, main="intrinsic MDS", ylim=c(-2,2), col=mylabs, pch=19)
#' plot(embed2ext, main="extrinsic MDS", ylim=c(-2,2), col=mylabs, pch=19)
#' par(opar)
#' 
#' @references 
#' \insertRef{torgerson_multidimensional_1952a}{Riemann}
#' 
#' @concept visualization
#' @export
riem.mds <- function(riemobj, ndim=2, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.mds : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim = max(2, round(ndim))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTE PAIRWISE DISTANCE
  distobj = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## COMPUTE MDS AND RETURN
  func.import = utils::getFromNamespace("hidden_cmds", "maotai")
  out.cmds    = func.import(distobj, ndim=myndim)
  return(out.cmds)
}