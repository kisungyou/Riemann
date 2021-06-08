#' Geodesic Interpolation of Multiple Points
#' 
#' Given 2 observations \eqn{X_1, X_2 \in \mathcal{M}}, find 
#' the interpolated points of a geodesic \eqn{\gamma(t)} for \eqn{t \in (0,1)} which 
#' assumes two endpoints \eqn{\gamma(0)=X_1} and \eqn{\gamma(1)=X_2}. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{2} manifold-valued data where the first object is the starting point.
#' @param vect a length-\eqn{T} vector in \eqn{(0,1)} for which the interpolations are taken.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return a 3d array where \eqn{T} slices along 3rd dimension are interpolated objects in matrix representation.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #       Geodesic Interpolation between (1,0) and (0,1) in S^1
#' #-------------------------------------------------------------------
#' ## PREPARE DATA
#' sp.start = c(1,0)
#' sp.end   = c(0,1)
#' sp.data  = wrap.sphere(rbind(sp.start, sp.end))
#' 
#' ## FIND THE INTERPOLATED POINT AT FOR t=0.1, 0.2, ..., 0.9.
#' myvect  = seq(from=0.1, to=0.9, by=0.1)
#' geo.int = riem.interps(sp.data, vect=myvect, geometry="intrinsic")
#' geo.ext = riem.interps(sp.data, vect=myvect, geometry="extrinsic")
#' 
#' geo.int = matrix(geo.int, byrow=TRUE, ncol=2) # re-arrange for plotting
#' geo.ext = matrix(geo.ext, byrow=TRUE, ncol=2)
#' 
#' ## VISUALIZE
#' #  Prepare Lines and Points
#' thetas  = seq(from=0, to=pi/2, length.out=100)
#' quarter = cbind(cos(thetas), sin(thetas))
#' 
#' pts.int = rbind(sp.start, geo.int, sp.end)
#' pts.ext = rbind(sp.start, geo.ext, sp.end)
#' col.int = c("black", rep("red",9),  "black")
#' col.ext = c("black", rep("blue",9), "black")
#' 
#' # Draw
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(quarter, main="intrinsic interpolation", # intrinsic geodesic
#'      xlab="x", ylab="y", type="l")
#' points(pts.int, col=col.int, pch=19)
#' for (i in 1:9){
#'   text(geo.int[i,1]*0.9, geo.int[i,2]*0.9, 
#'        paste0(round(i/10,2)), col="red")
#' }
#' plot(quarter, main="extrinsic interpolation", # intrinsic geodesic
#'      xlab="x", ylab="y", type="l")
#' points(pts.ext, col=col.ext, pch=19)
#' for (i in 1:9){
#'   text(geo.ext[i,1]*0.9, geo.ext[i,2]*0.9, 
#'        paste0(round(i/10,2)), col="blue")
#' }
#' par(opar)
#' 
#' @concept basic
#' @export
riem.interps <- function(riemobj, vect=c(0.25, 0.5, 0.75), geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.interps : input ",DNAME," should be an object of 'riemdata' class."))
  }
  if (length(riemobj$data)!=2){
    stop("* riem.interps : for 'riem.interps', please have two objects only.")
  }
  if ((!is.vector(vect))||(any(vect<=0))||(any(vect>=1))){
    stop("* riem.interps : 'vect' should contain numbers in (0,1).")
  }
  if (length(vect)==1){
    stop("* riem.interps : 'vect' should contain multiple 't' values. For a single number, use 'riem.interp' instead.")
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTATION
  outmat = basic_interpolate(riemobj$name, mygeom, riemobj$data[[1]], riemobj$data[[2]], vect)
  return(outmat)
}
