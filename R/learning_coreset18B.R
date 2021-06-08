#' Build Lightweight Coreset
#' 
#' Given manifold-valued data \eqn{X_1,X_2,\ldots,X_N \in \mathcal{M}}, this algorithm 
#' finds the coreset of size \eqn{M} that can be considered as a compressed representation 
#' according to the lightweight coreset construction scheme proposed by the reference below.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param M the size of coreset (default: \eqn{N/2}).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{coreid}{a length-\eqn{M} index vector of the coreset.}
#' \item{weight}{a length-\eqn{M} vector of weights for each element.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # * 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # * 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # * 10 perturbed data points near (0,0,1) on S^2 in R^3
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
#' 
#' ## MDS FOR VISUALIZATION
#' embed2 = riem.mds(myriem, ndim=2)$embed
#' 
#' ## FIND CORESET OF SIZES 3, 6, 9
#' core1 = riem.coreset18B(myriem, M=3)
#' core2 = riem.coreset18B(myriem, M=6)
#' core3 = riem.coreset18B(myriem, M=9)
#' 
#' col1 = rep(1,30); col1[core1$coreid] = 2
#' col2 = rep(1,30); col2[core2$coreid] = 2
#' col3 = rep(1,30); col3[core3$coreid] = 2
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embed2, pch=19, col=col1, main="coreset size=3")
#' plot(embed2, pch=19, col=col2, main="coreset size=6")
#' plot(embed2, pch=19, col=col3, main="coreset size=9")
#' par(opar)
#' 
#' @references 
#' \insertRef{bachem_scalable_2018a}{Riemann}
#' 
#' @concept learning
#' @export
riem.coreset18B <- function(riemobj, M=length(riemobj$data)/2, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE : EXPLICIT
  N     = length(riemobj$data)
  mygeo = ifelse(missing(geometry),"intrinsic",match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  mym   = max(2, round(M))
  
  ## PREPARE : IMPLICIT
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 50))
  myeps  = min(1e-5, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-5)))
  
  ## RUN
  cpprun = learning_coreset18B(riemobj$name, mygeo, riemobj$data, mym, myiter, myeps)
  coreids = as.vector(cpprun$id)+1
  probvec = as.vector(cpprun$qx)
    
  ## WRAP AND RETURN
  output = list()
  output$coreid = coreids
  output$weight = 1/(probvec[coreids]*mym)
  return(output)
}