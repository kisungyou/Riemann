#' K-Means Clustering
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform k-means clustering by minimizing within-cluster sum of squares (WCSS). 
#' Since the problem is NP-hard and sensitive to the initialization, we provide an 
#' option with multiple starts and return the best result with respect to WCSS.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{algorithm}{(case-insensitive) name of an algorithm; \code{"MacQueen"} (default), or \code{"Lloyd"}.}
#' \item{init}{(case-insensitive) name of an initialization scheme; \code{"plus"} for k-means++ (default), or \code{"random"}.}
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
#' ## K-MEANS WITH K=2,3,4
#' clust2 = riem.kmeans(myriem, k=2)
#' clust3 = riem.kmeans(myriem, k=3)
#' clust4 = riem.kmeans(myriem, k=4)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(mds2d, pch=19, main="true label", col=mylabs)
#' plot(mds2d, pch=19, main="K=2", col=clust2$cluster)
#' plot(mds2d, pch=19, main="K=3", col=clust3$cluster)
#' plot(mds2d, pch=19, main="K=4", col=clust4$cluster)
#' par(opar)
#' 
#' @seealso \code{\link{riem.kmeanspp}}
#' 
#' @references 
#' \insertRef{lloyd_least_1982}{Riemann}
#' 
#' \insertRef{macqueen_methods_1967}{Riemann}
#' 
#' @concept clustering
#' @export
riem.kmeans <- function(riemobj, k=2, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE
  N          = length(riemobj$data)
  par.geo    = ifelse(missing(geometry),"intrinsic",match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  par.k      = max(1, round(k))
  
  # IMPLICIT PARAMETERS 
  pars   = list(...)
  pnames = names(pars)
  par.iter = ifelse(("maxiter"%in%pnames), max(50, round(pars$maxiter)), 50)
  par.init = ifelse(("init"%in%pnames), match.arg(tolower(pars$init),c("plus","random")), "plus")
  par.alg  = ifelse(("algorithm"%in%pnames), match.arg(tolower(pars$algorithm),c("macqueen","lloyd")), "macqueen")
  par.nstart = ifelse(("nstart"%in%pnames), max(2, round(pars$nstart)), 5)

  ## INITIALIZATION
  rec.lab0 = list()
  if (all(par.init=="random")){
    for (i in 1:par.nstart){
      rec.lab0[[i]] = base::sample(c(c(1:par.k), sample(1:par.k, (N-par.k), replace = TRUE)))
    }  
  } else {
    distobj     = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, par.geo))
    func.import = utils::getFromNamespace("hidden_kmeanspp", "maotai")
    for (i in 1:par.nstart){
      rec.lab0[[i]] = as.vector(as.integer(func.import(distobj, k=par.k)$cluster))
    }
  }

  ## MAIN RUN : RETURN THE BEST WITH RESPECT TO WCSS
  rec.list = list()
  rec.SSE  = rep(0,par.nstart)
  
  for (i in 1:par.nstart){
    rec.list[[i]] = switch(par.alg,
                           macqueen = clustering_kmeans_macqueen(riemobj$name, par.geo, riemobj$data, par.iter, 0.01, rec.lab0[[i]]),
                           lloyd    = clustering_kmeans_lloyd(riemobj$name, par.geo, riemobj$data, par.iter, 0.01, rec.lab0[[i]]))
    rec.SSE[i] = rec.list[[i]]$WCSS
  }
  
  ## SELECT, WRAP, AND RETURN
  bestout = rec.list[[which.min(rec.SSE)]]
  
  output = list()
  output$cluster = as.vector(bestout$label)+1
  output$means   = bestout$means
  output$score   = as.double(bestout$WCSS)
  return(output)
}