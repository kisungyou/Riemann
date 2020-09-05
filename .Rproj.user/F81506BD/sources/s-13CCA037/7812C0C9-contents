#' Compute Pairwise Distances for Two Sets of Data
#' 
#' Given \eqn{M} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}} and 
#' \eqn{N} observations \eqn{Y_1, Y_2, \ldots, Y_N \in \mathcal{M}}, 
#' compute pairwise distances between two sets' elements.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return an \eqn{(M\times N)} matrix of distances.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with two types
#' #
#' #  group1 : 10 perturbed data points near (0,0,1) on S^2 in R^3
#' #  group2 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#' #           10 perturbed data points near (0,1,0) on S^2 in R^3
#' #           10 perturbed data points near (0,0,1) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata1 = list()
#' mydata2 = list()
#' for (i in 1:10){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata1[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 1:10){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 21:30){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem1 = wrap.sphere(mydata1)
#' myriem2 = wrap.sphere(mydata2)
#' 
#' ## COMPARE TWO DISTANCES
#' dint = riem.pdist2(myriem1, myriem2, geometry="intrinsic")
#' dext = riem.pdist2(myriem1, myriem2, geometry="extrinsic")
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' image(dint[nrow(dint):1,], main="intrinsic", axes=FALSE)
#' image(dext[nrow(dext):1,], main="extrinsic", axes=FALSE)
#' par(opar)
#' 
#' @concept basic
#' @export
riem.pdist2 <- function(riemobj1, riemobj2, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* riem.pdist2 : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* riem.pdist2 : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  if (riemobj1$name!=riemobj2$name){
    stop("* riem.pdist2 : two inputs are from different manifolds.")
  }
  if (!all(riemobj1$size==riemobj2$size)){
    stop("* riem.pdist2 : two inputs are of different size.")
  }
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))

  ## COMPUTE & RETURN
  distmat = basic_pdist2(riemobj1$name, riemobj1$data, riemobj2$data, mygeometry)
  return(distmat)
}
