#' K-Medoids Clustering
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform k-medoids clustering using pairwise distances.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return a named list containing\describe{
#' \item{medoids}{a length-\eqn{k} vector of medoids' indices.}
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
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
#' ## K-MEDOIDS WITH K=2,3,4
#' clust2 = riem.kmedoids(myriem, k=2)
#' clust3 = riem.kmedoids(myriem, k=3)
#' clust4 = riem.kmedoids(myriem, k=4)
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
#' @seealso \code{\link[cluster]{pam}} 
#' @concept clustering
#' @export
riem.kmedoids <- function(riemobj, k=2, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.kmedoids : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myk    = max(0, round(k))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTE PAIRWISE DISTANCE
  distobj = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## RUN K-MEDOIDS
  func.import  = utils::getFromNamespace("hidden_kmedoids", "maotai")
  obj.kmedoids = func.import(distobj, nclust=myk) 
  
  ## WRAP AND RETURN
  output = list()
  output$medoids = obj.kmedoids$id.med
  output$cluster = as.integer(obj.kmedoids$clustering)
  return(obj.kmedoids)
  
}