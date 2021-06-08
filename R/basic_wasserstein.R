#' Wasserstein Distance between Empirical Measures
#' 
#' Given two empirical measures \eqn{\mu, \nu} consisting of \eqn{M} and \eqn{N} observations, \eqn{p}-Wasserstein distance for \eqn{p\geq 1} between two empirical measures 
#' is defined as 
#' \deqn{\mathcal{W}_p (\mu, \nu) = \left( \inf_{\gamma \in \Gamma(\mu, \nu)} \int_{\mathcal{M}\times \mathcal{M}} d(x,y)^p d \gamma(x,y) \right)^{1/p}}
#' where \eqn{\Gamma(\mu, \nu)} denotes the collection of all measures/couplings on \eqn{\mathcal{M}\times \mathcal{M}} 
#' whose marginals are \eqn{\mu} and \eqn{\nu} on the first and second factors, respectively.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data, which are atoms of \eqn{\mu}.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data, which are atoms of \eqn{\nu}.
#' @param p an exponent for Wasserstein distance \eqn{\mathcal{W}_p} (default: 2).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{weight1}{a length-\eqn{M} weight vector for \eqn{\mu}; if \code{NULL} (default), uniform weight is set.}
#' \item{weight2}{a length-\eqn{N} weight vector for \eqn{\nu}; if \code{NULL} (default), uniform weight is set.}
#' }
#' 
#' @return a named list containing \describe{
#' \item{distance}{\eqn{\mathcal{W_p}} distance between two empirical measures.}
#' \item{plan}{an \eqn{(M\times N)} matrix whose rowSums and columnSums are \code{weight1} and \code{weight2} respectively.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with two types
#' #
#' # class 1 : 20 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 30 perturbed data points near (0,1,0) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata1 = list()
#' mydata2 = list()
#' for (i in 1:20){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata1[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 1:30){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem1 = wrap.sphere(mydata1)
#' myriem2 = wrap.sphere(mydata2)
#' 
#' ## COMPUTE p-WASSERSTEIN DISTANCES
#' dist1 = riem.wasserstein(myriem1, myriem2, p=1)
#' dist2 = riem.wasserstein(myriem1, myriem2, p=2)
#' dist5 = riem.wasserstein(myriem1, myriem2, p=5)
#' 
#' pm1 = paste0("p=1: dist=",round(dist1$distance,3))
#' pm2 = paste0("p=2: dist=",round(dist2$distance,3))
#' pm5 = paste0("p=5: dist=",round(dist5$distance,3))
#' 
#' ## VISUALIZE TRANSPORT PLAN AND DISTANCE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(dist1$plan, axes=FALSE, main=pm1)
#' image(dist2$plan, axes=FALSE, main=pm2)
#' image(dist5$plan, axes=FALSE, main=pm5)
#' par(opar)
#' 
#' @concept basic
#' @export
riem.wasserstein <- function(riemobj1, riemobj2, p=2, geometry=c("intrinsic","extrinsic"), ...){
  ## INPUTS : EXPLICIT
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* riem.wasserstein : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* riem.wasserstein : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  myp = max(1, as.double(p))
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  M = length(riemobj1$data)
  N = length(riemobj2$data)
  
  ## INPUTS : IMPLICIT
  param  = list(...)
  pnames = names(param)
  
  if ("weight1"%in%pnames){
    myweight1 = param$weight1
    if ((length(myweight1)<1)&&is.null(myweight1)){
      myweight1 = rep(1/M, M)
    } else {
      myweight1 = myweight1/sum(myweight1)
    }
  } else {
    myweight1 = rep(1/M, M)
  }
  if ("weight2"%in%pnames){
    myweight2 = param$weight2
    if ((length(myweight2)<1)&&is.null(myweight2)){
      myweight2 = rep(1/N, N)
    } else {
      myweight2 = myweight2/sum(myweight2)
    }
  } else {
    myweight2 = rep(1/N, N)
  }
  
  if ((length(myweight1)!=M)||(any(myweight1<0))){
    stop("* riem.wasserstein : 'weight1' should be of length matching to 'riemobj1' & no negative values are admitted.")
  }
  if ((length(myweight2)!=N)||(any(myweight2<0))){
    stop("* riem.wasserstein : 'weight2' should be of length matching to 'riemobj2' & no negative values are admitted.")
  }
  
  ## SWITCHING, COMPUTATION, AND RETURN
  dxy    = riem.pdist2(riemobj1, riemobj2, geometry=mygeometry)
  # output = T4transport::wassersteinD(dxy, myp, wx=myweight1, wy=myweight2)
  output = T4transport_wassersteinD(dxy, myp, wx=myweight1, wy=myweight2)
  return(output)
}

# ## WANT TO SEE CONCENTRATION OF EMPIRICAL DISTANCE FOR
# ## TWO EMPIRICAL MEASURES : EXPECTED DISTANCE IS "pi/2"
# 
# niter = 10000
# distn = rep(0,niter)
# for (it in 1:niter){
#   mydata1 = list()
#   mydata2 = list()
#   for (i in 1:30){
#     tgt = c(1, stats::rnorm(2, sd=0.1))
#     mydata1[[i]] = tgt/sqrt(sum(tgt^2))
#   }
#   for (i in 1:30){
#     tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#     mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#   }
#   myriem1  = wrap.sphere(mydata1)
#   myriem2  = wrap.sphere(mydata2)
#   distn[it] = riem.wasserstein(myriem1, myriem2)$distance
#   if (it%%100 == 0){
#     print(paste0("iteration ",it,"/",niter," complete.."))
#   }
# }
# 
# # Visualize
# opar <- par(no.readonly=TRUE)
# hist(distn, main="Monte Carlo Distribution")
# abline(v=pi/2, lwd=2, col="red")
# par(opar)