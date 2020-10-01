#' Sammon Mapping
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' apply Sammon mapping, a non-linear dimensionality reduction method. Since 
#' the method depends only on the pairwise distances of the data, it can be 
#' adapted to the manifold-valued data.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
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
#' ## COMPARE SAMMON WITH MDS
#' embed2mds = riem.mds(myriem, ndim=2)$embed
#' embed2sam = riem.sammon(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2mds, col=mylabs, pch=19, main="MDS")
#' plot(embed2sam, col=mylabs, pch=19, main="Sammon mapping")
#' par(opar)
#' 
#' @references 
#' \insertRef{sammon_nonlinear_1969a}{Riemann}
#' 
#' @concept visualization
#' @export
riem.sammon <- function(riemobj, ndim=2, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.sammon : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim = max(1, round(ndim))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  # IMPLICIT PARAMETERS 
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 50))
  myeps  = min(1e-5, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-5)))
  
  ## RUN FROM RCPP
  return(visualize_sammon(riemobj$name, mygeom, riemobj$data, myndim, myiter, myeps))
}