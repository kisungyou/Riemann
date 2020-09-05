#' Find K-Nearest Neighbors
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' \code{riem.knn} constructs \eqn{k}-nearest neighbors. This is 
#' a wrapper for a main function in \pkg{nabor} package. Note that the original 
#' function contains index for each data point itself, but our function does not 
#' consider self-neighborhood scenario. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of neighbors to find.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters to be passed onto \code{\link[nabor]{knn}} function.
#' 
#' @return a named list containing\describe{
#' \item{nn.idx}{an \eqn{(N \times k)} neighborhood index matrix.}
#' \item{nn.dists}{an \eqn{(N\times k)} distances from a point to its neighbors.}
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
#' mylabs = rep(c(2,3,4), each=10)
#' 
#' ## K-NN CONSTRUCTION WITH K=5 & K=10
#' knn1 = riem.knn(myriem, k=5)
#' knn2 = riem.knn(myriem, k=10)
#' 
#' ## MDS FOR VISUALIZATION
#' embed2 = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2, pch=19, main="knn with k=4", col=mylabs)
#' for (i in 1:30){
#'   for (j in 1:5){
#'     lines(embed2[c(i,knn1$nn.idx[i,j]),])
#'   }
#' }
#' plot(embed2, pch=19, main="knn with k=8", col=mylabs)
#' for (i in 1:30){
#'   for (j in 1:10){
#'     lines(embed2[c(i,knn2$nn.idx[i,j]),])
#'   }
#' }
#' par(opar)
#' 
#' @seealso \code{\link[nabor]{knn}}
#' @concept learning
#' @export
riem.knn <- function(riemobj, k=2, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.knn : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myk    = max(0, round(k))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTE PAIRWISE DISTANCE
  distobj = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## RUN KNN
  func.import = utils::getFromNamespace("hidden_knn", "maotai")
  obj.knn     = func.import(distobj, nnbd=(myk+1), ...)
  
  ## WRAP AND RETURN
  output = list()
  output$nn.idx   = obj.knn$nn.idx[,2:(myk+1)]
  output$nn.dists = obj.knn$nn.dists[,2:(myk+1)]
  return(output)
}