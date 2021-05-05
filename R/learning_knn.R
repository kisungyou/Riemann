#' Find K-Nearest Neighbors
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' \code{riem.knn} constructs \eqn{k}-nearest neighbors.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of neighbors to find.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
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
#' @concept learning
#' @export
riem.knn <- function(riemobj, k=2, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.knn : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myk    = max(0, round(k))
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTE PAIRWISE DISTANCE
  distobj = as.matrix(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## COMPUTE AND RETURN
  return(nearest_neighbor(distobj, myk))
}


#' @keywords internal
#' @noRd
nearest_neighbor <- function(dmat, k){
  n = base::nrow(dmat)
  
  nn.idx   = array(0,c(n,k))
  nn.dists = array(0,c(n,k))
  
  for (i in 1:n){
    tgt  = as.vector(dmat[i,])
    i_id = order(tgt)[2:(k+1)]

    nn.idx[i,]   = i_id
    nn.dists[i,] = tgt[i_id]
  }
  
  output = list()
  output$nn.idx   = nn.idx
  output$nn.dists = nn.dists
  return(output)
}


# library(usmap)
# library(ggplot2)
# data("cities")
# mygeo = usmap_transform(data.frame(lon=cities$coord[,2], lat=cities$coord[,1]))
# 
# myriem = riem.sphere(cities$cartesian)
# 
# plot_usmap(regions="states") + 
#   geom_point(data=mygeo, aes(x=lon.1, y=lat.1), alpha=0.25)
