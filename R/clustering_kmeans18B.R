#' K-Means Clustering with Lightweight Coreset
#' 
#' The modified version of lightweight coreset for scalable \eqn{k}-means computation 
#' is applied for manifold-valued data \eqn{X_1,X_2,\ldots,X_N \in \mathcal{M}}. 
#' The smaller the set is, the faster the execution becomes with potentially larger quantization errors.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param M the size of coreset (default: \eqn{N/2}).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{nstart}{the number of random starts (default: 5).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
#' \item{means}{a 3d array where each slice along 3rd dimension is a matrix representation of class mean.}
#' \item{score}{within-cluster sum of squares (WCSS).}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # class 1 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # class 3 : 10 perturbed data points near (0,0,1) on S^2 in R^3
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
#' ## TRY DIFFERENT SIZES OF CORESET WITH K=4 FIXED
#' core1 = riem.kmeans18B(myriem, k=3, M=5)
#' core2 = riem.kmeans18B(myriem, k=3, M=10)
#' core3 = riem.kmeans18B(myriem, k=3, M=15)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(mds2d, pch=19, main="true label", col=mylabs)
#' plot(mds2d, pch=19, main="kmeans18B: M=5",  col=core1$cluster)
#' plot(mds2d, pch=19, main="kmeans18B: M=10", col=core2$cluster)
#' plot(mds2d, pch=19, main="kmeans18B: M=15", col=core3$cluster)
#' par(opar)
#' 
#' @references 
#' \insertRef{bachem_scalable_2018a}{Riemann}
#' 
#' @seealso \code{\link{riem.coreset18B}}
#' @concept clustering
#' @export
riem.kmeans18B <- function(riemobj, k=2, M=length(riemobj$data)/2, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE
  N       = length(riemobj$data)
  par.geo = ifelse(missing(geometry),"intrinsic",match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  par.k   = max(1, round(k))
  par.m   = max(2*par.k, round(M))
  
  # IMPLICIT PARAMETERS 
  pars   = list(...)
  pnames = names(pars)
  par.iter   = ifelse(("maxiter"%in%pnames), max(50, round(pars$maxiter)), 50)
  par.nstart = ifelse(("nstart"%in%pnames), round(pars$nstart), 5)
  
  ## MAIN RUN : MULTIPLE START WITH BEST WCSS
  rec.list = list()
  rec.SSE  = rep(0,par.nstart)
  
  for (i in 1:par.nstart){
    rec.list[[i]] = clustering_kmeans18B(riemobj$name, par.geo, riemobj$data, par.k, par.m, par.iter)
    rec.SSE[i]    = rec.list[[i]]$wcss
  }
  
  ## SELECT, WRAP, AND RETURN
  bestout = rec.list[[which.min(rec.SSE)]]
  output  = list()
  output$cluster = as.vector(bestout$cluster)+1
  output$means   = bestout$means
  output$score   = as.double(bestout$wcss)
  return(output)
}