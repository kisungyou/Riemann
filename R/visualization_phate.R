#' PHATE
#' 
#' PHATE is a nonlinear manifold learning method that is specifically targeted at 
#' improving diffusion maps by incorporating data-adaptive kernel construction, 
#' detection of optimal time scale, and information-theoretic metric measures.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters for \code{PHATE} including \describe{
#' \item{nbdk}{size of nearest neighborhood (default: 5).}
#' \item{alpha}{decay parameter for Gaussian kernel exponent (default: 2).}
#' \item{potential}{type of potential distance transformation; \code{"log"} or \code{"sqrt"} (default: \code{"log"}).}
#' }
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' }
#' 
#' @examples 
#' \donttest{
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
#' ## PHATE EMBEDDING WITH LOG & SQRT POTENTIAL 
#' phate_log  = riem.phate(myriem, potential="log")$embed
#' phate_sqrt = riem.phate(myriem, potential="sqrt")$embed
#' embed_mds  = riem.mds(myriem)$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embed_mds,  col=mylabs, pch=19, main="MDS" )
#' plot(phate_log,  col=mylabs, pch=19, main="PHATE+Log")
#' plot(phate_sqrt, col=mylabs, pch=19, main="PHATE+Sqrt")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{moon_visualizing_2019}{Riemann}
#' 
#' @concept visualization
#' @export
riem.phate <- function(riemobj, ndim=2, geometry=c("intrinsic","extrinsic"), ...){
  ## INPUT : EXPLICIT
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.phate : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim = max(2, round(ndim))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("nbdk"%in%pnames){
    myk = max(1, round(params$nbdk))
  } else {
    myk = 5
  }
  if ("alpha"%in%pnames){
    myalpha = max(.Machine$double.eps, as.double(params$alpha))
  } else {
    myalpha = 2
  }
  if ("potential"%in%pnames){
    mydtype = match.arg(params$potential, c("log","sqrt"))
  } else {
    mydtype = "log"
  }
  
  ## COMPUTE PAIRWISE DISTANCES
  distobj = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## POTENTIAL
  phate_op  <- utils::getFromNamespace("hidden_PHATE","maotai")
  phate_run <- phate_op(distobj, nbdk=myk, alpha=myalpha)
  phate_t   <- phate_run$t
  
  ## TRANSFORM & MMDS
  if (all(mydtype=="sqrt")){
    optP = base::sqrt(phate_run$P)
  } else {
    optP = base::log(phate_run$P + (1e-8))
  }
  optPdist = stats::as.dist(cpp_pdist(optP))
  mmds_op  = utils::getFromNamespace("hidden_mmds","maotai")
  
  ## RETURN
  output = list()
  output$embed = mmds_op(optPdist, ndim=myndim, abstol=1e-8)
  return(output)
}


