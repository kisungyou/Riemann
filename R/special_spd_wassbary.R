#' Wasserstein Barycenter of SPD Matrices
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N} in SPD manifold, compute 
#' the \eqn{L_2}-Wasserstein barycenter that minimizes
#' \deqn{\sum_{n=1}^N \lambda_i \mathcal{W}_2 (N(X), N(X_i))^2}
#' where \eqn{N(X)} denotes the zero-mean Gaussian measure with covariance \eqn{X}.
#' 
#' @param spdobj a S3 \code{"riemdata"} class of SPD-valued data of \eqn{(p\times p)} matrices.
#' @param weight weight of observations; if \code{NULL} it assumes equal weights, or a nonnegative length-\eqn{N} vector that sums to 1 should be given.
#' @param method name of the algortihm to be used; one of the \code{"RU02"}, \code{"AE16"}.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:20).}
#' \item{abstol}{tolerance level for stopping criterion (default: 1e-8).}
#' }
#' 
#' @return a \eqn{(p\times p)} Wasserstein barycenter matrix.
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #        Covariances from standard multivariate Gaussians.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' ndata = 20
#' pdim  = 10
#' mydat = array(0,c(pdim,pdim,ndata))
#' for (i in 1:ndata){
#'   mydat[,,i] = stats::cov(matrix(rnorm(100*pdim), ncol=pdim))
#' }
#' myriem = wrap.spd(mydat)
#' 
#' ## COMPUTE BY DIFFERENT ALGORITHMS
#' baryRU <- spd.wassbary(myriem, method="RU02")
#' baryAE <- spd.wassbary(myriem, method="AE16")
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(diag(pdim), axes=FALSE, main="True Covariance")
#' image(baryRU, axes=FALSE, main="by RU02")
#' image(baryAE, axes=FALSE, main="by AE16")
#' par(opar)
#' }
#' 
#' 
#' @concept spd
#' @export
spd.wassbary <- function(spdobj, weight=NULL, method=c("RU02","AE16"), ...){
  # INPUTS : EXPLICIT
  DNAME = paste0("'",deparse(substitute(spdobj)),"'")
  if ((!inherits(spdobj,"riemdata"))||(!all(spdobj$name=="spd"))){
    stop(paste0("* spd.wassbary : input ",DNAME," should be an object of 'riemdata' class on 'spd' manifold.."))
  }
  N = length(spdobj$data)
  if ((length(weight)<1)&&(is.null(weight))){
    myweight = rep(1/N, N)
  } else {
    myweight = check_weight(weight, N, "spd.wassbary")
  }
  mymethod = match.arg(method)
  
  # INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  if ("maxiter"%in%pnames){
    maxiter = max(5, round(params$maxiter))
  } else {
    maxiter = 20
  }
  if ("abstol"%in%pnames){
    abstol = max(100*.Machine$double.eps, as.double(params$abstol))
  } else {
    abstol = 10^(-8)
  }

  # COMPUTE
  output = switch(mymethod,
                  "RU02" = spdwass_baryRU02(spdobj$data, myweight, maxiter, abstol),
                  "AE16" = spdwass_baryAE16(spdobj$data, myweight, maxiter, abstol))
  return(output)
}
  