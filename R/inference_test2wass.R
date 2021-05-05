#' Two-Sample Test with Wasserstein Metric
#' 
#' Given \eqn{M} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}} and 
#' \eqn{N} observations \eqn{Y_1, Y_2, \ldots, Y_N \in \mathcal{M}}, permutation 
#' test based on the Wasserstein metric (see \code{\link{riem.wasserstein}} for 
#' more details) is applied to test whether two distributions are same or not, i.e.,
#' \deqn{H_0~:~\mathcal{P}_X = \mathcal{P}_Y}
#' with Wasserstein metric \eqn{\mathcal{W}_p} being the measure of discrepancy 
#' between two samples.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param p an exponent for Wasserstein distance \eqn{\mathcal{W}_p} (default: 2).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{nperm}{the number of permutations (default: 999).}
#' \item{use.smooth}{a logical; \code{TRUE} to use a smoothed Wasserstein distance, \code{FALSE} otherwise.}
#' }
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
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
#' for (i in 1:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem1 = wrap.sphere(mydata1)
#' myriem2 = wrap.sphere(mydata2)
#' 
#' ## PERFORM PERMUTATION TEST
#' #  it is expected to return a very small number, but 
#' #  small number of 'nperm' may not give a reasonable p-value.
#' \donttest{
#' riem.test2wass(myriem1, myriem2, nperm=99, use.smooth=FALSE)
#' }
#' 
#' \dontrun{
#' ## CHECK WITH EMPIRICAL TYPE-1 ERROR
#' set.seed(777)
#' ntest = 1000
#' pvals = rep(0,ntest)
#' 
#' for (i in 1:ntest){
#'   X = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Y = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Xnorm = X/sqrt(rowSums(X^2))
#'   Ynorm = Y/sqrt(rowSums(Y^2))
#'   
#'   Xriem = wrap.sphere(Xnorm)
#'   Yriem = wrap.sphere(Ynorm)
#'   pvals[i] = riem.test2wass(Xriem, Yriem, nperm=999)$p.value
#'   print(paste0("iteration ",i,"/",ntest," complete.."))
#' }
#' 
#' emperr = round(sum((pvals <= 0.05))/ntest, 5)
#' print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
#' }
#' 
#' @concept inference
#' @export
riem.test2wass <- function(riemobj1, riemobj2, p=2, geometry=c("intrinsic","extrinsic"), ...){
  ## INPUTS : EXPLICIT
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* wass.test2wass : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* wass.test2wass : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  myp = max(1, as.double(p))
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  M = length(riemobj1$data)
  N = length(riemobj2$data)
  w1 = rep(1/M, M)
  w2 = rep(1/N, N)
  
  ## INPUTS : IMPLICIT
  param    = list(...)
  pnames   = names(param)
  mynperm  = ifelse(("nperm"%in%pnames), max(9, round(param$nperm)), 999)
  myipot   = as.logical(ifelse(("use.smooth"%in%pnames), param$use.smooth, TRUE))
  
  ## COMPUTE : DISTANCE AND STATISTIC UNDER NULL
  distmat = basic_pdist(riemobj1$name, c(riemobj1$data, riemobj2$data), mygeometry)
  if (myipot){
    thestat = T4transport_ipotD(distmat[1:M, (M+1):(M+N)],p=myp,wx = w1, wy=w2)$distance
  } else {
    thestat = T4transport_wassersteinD(distmat[1:M, (M+1):(M+N)], myp, wx=w1, wy=w2)$distance
  }
  
  ## COMPUTE : ITERATION
  distvals = rep(0, mynperm)
  for (i in 1:mynperm){
    id.all  = sample(1:(M+N))
    id.gp1  = id.all[1:M]
    id.gp2  = id.all[(M+1):(M+N)]
    partdxy = distmat[id.gp1, id.gp2]
    
    if (myipot){
      distvals[i] = T4transport_ipotD(partdxy,p=myp,wx = w1, wy=w2)$distance
    } else {
      distvals[i] = T4transport_wassersteinD(partdxy, myp, wx=w1, wy=w2)$distance
    }
  }
  
  ## WRAP
  pvalue   = (sum(distvals >= thestat)+1)/(mynperm+1)
  dataname = paste0(DNAME1," and ",DNAME2)
  mfdname  = wrap_mfd2full(riemobj1$name)
  hname    = paste0("Wasserstein Two-Sample Test on ",mfdname," Manifold")
  Ha       = "two distributions are not equal."
  names(thestat) = "Wmn"
  
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name=dataname)
  class(res) = "htest"
  return(res)
}