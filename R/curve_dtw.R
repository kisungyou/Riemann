#' Dynamic Time Warping Distance 
#' 
#' Given two time series - a query \eqn{X = (X_1,X_2,\ldots,X_N)} and a reference \eqn{Y = (Y_1,Y_2,\ldots,Y_M)}, 
#' \code{riem.dtw} computes the most basic version of Dynamic Time Warping (DTW) distance between two series using a symmetric step pattern, meaning 
#' no window constraints and others at all. Although the scope of DTW in Euclidean space-valued objects is rich, it is scarce for manifold-valued curves. 
#' If you are interested in the topic, we refer to \pkg{dtw} package.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data along the curve.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data along the curve.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return the distance value.
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #                          Curves on Sphere
#' #
#' #  curve1 : y = 0.5*cos(x) on the tangent space at (0,0,1)
#' #  curve2 : y = 0.5*sin(x) on the tangent space at (0,0,1)
#' # 
#' #  we will generate two sets for curves of different sizes.
#' #-------------------------------------------------------------------
#' ## GENERATION
#' clist = list()
#' for (i in 1:10){ # curve type 1
#'   vecx = seq(from=-0.9, to=0.9, length.out=sample(10:50, 1))
#'   vecy = 0.5*cos(vecx) + rnorm(length(vecx), sd=0.1)
#'   mats = cbind(vecx, vecy, 1)
#'   clist[[i]] = wrap.sphere(mats/sqrt(rowSums(mats^2)))
#' }
#' for (i in 1:10){ # curve type 2
#'   vecx = seq(from=-0.9, to=0.9, length.out=sample(10:50, 1))
#'   vecy = 0.5*sin(vecx) + rnorm(length(vecx), sd=0.1)
#'   mats = cbind(vecx, vecy, 1)
#'   clist[[i+10]] = wrap.sphere(mats/sqrt(rowSums(mats^2)))
#' }
#' 
#' ## COMPUTE DISTANCES
#' outint = array(0,c(20,20))
#' outext = array(0,c(20,20))
#' for (i in 1:19){
#'   for (j in 2:20){
#'     outint[i,j] <- outint[j,i] <- riem.dtw(clist[[i]], clist[[j]], 
#'                                            geometry="intrinsic")
#'     outext[i,j] <- outext[j,i] <- riem.dtw(clist[[i]], clist[[j]],
#'                                            geometry="extrinsic")
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(outint[,20:1], axes=FALSE, main="intrinsic DTW Distance")
#' image(outext[,20:1], axes=FALSE, main="extrinsic DTW Distance")
#' par(opar)
#' }
#' 
#' @concept curve
#' @export
riem.dtw <- function(riemobj1, riemobj2, geometry=c("intrinsic","extrinsic")){
  ## PREPARE : EXPLICIT
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* riem.dtw : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* riem.dtw : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  if (!all(riemobj1$name==riemobj2$name)){
    stop("* riem.dtw : two inputs are from different manifolds.")
  }
  if (!all(riemobj1$size==riemobj2$size)){
    stop("* riem.dtw : two inputs are of different size.")
  }
  mymfd   = riemobj1$name
  mygeom  = ifelse(missing(geometry),"intrinsic",
                   match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  # COMPUTE THE DISTANCE
  outdist = curvedist_dtwbasic(mymfd, mygeom, riemobj1$data, riemobj2$data) # dtw:: step pattern of symmetric1
  return(outdist)
}