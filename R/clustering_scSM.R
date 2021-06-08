#' Spectral Clustering by Shi and Malik (2000)
#' 
#' The version of Shi and Malik first constructs the affinity matrix
#' \deqn{A_{ij} = \exp(-d(x_i, d_j)^2 / \sigma^2)}
#' where \eqn{\sigma} is a common bandwidth parameter and performs k-means clustering on 
#' the row-space of eigenvectors for the random-walk graph laplacian matrix
#' \deqn{L=D^{-1}(D-A)}.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters (default: 2).
#' @param sigma bandwidth parameter (default: 1).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return a named list containing 
#' \describe{
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).} 
#' \item{eigval}{eigenvalues of the graph laplacian's spectral decomposition.}
#' \item{embeds}{an \eqn{(N\times k)} low-dimensional embedding.}
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
#' lab    = rep(c(1,2,3), each=10)
#' 
#' ## CLUSTERING WITH DIFFERENT K VALUES
#' cl2 = riem.scSM(myriem, k=2)$cluster
#' cl3 = riem.scSM(myriem, k=3)$cluster
#' cl4 = riem.scSM(myriem, k=4)$cluster
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(mds2d, col=lab, pch=19, main="true label")
#' plot(mds2d, col=cl2, pch=19, main="riem.scSM: k=2")
#' plot(mds2d, col=cl3, pch=19, main="riem.scSM: k=3")
#' plot(mds2d, col=cl4, pch=19, main="riem.scSM: k=4")
#' par(opar)
#' 
#' @references 
#' Shi J, Malik J (2000). “Normalized Cuts and Image Segmentation." \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 22(8):888–905.
#' 
#' @concept clustering
#' @export
riem.scSM <- function(riemobj, k=2, sigma=1, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.scSM : input ",DNAME," should be an object of 'riemdata' class."))
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  myk   = max(1, round(k))
  mysig = max(sqrt(.Machine$double.eps), as.double(sigma))
  
  ## COMPUTE DISTANCE
  pdmat   = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  
  ## RUN SPECTRAL CLUSTERING
  runT4cluster = T4cluster::scSM(pdmat, k=myk, sigma=mysig)

  ## WRAP AND RETURN
  output = list()
  output$cluster = runT4cluster$cluster
  output$eigval  = runT4cluster$eigval
  output$embeds  = runT4cluster$embeds
  return(output)
}