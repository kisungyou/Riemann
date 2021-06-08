#' Geodesic Interpolation
#' 
#' Given 2 observations \eqn{X_1, X_2 \in \mathcal{M}}, find the interpolated 
#' point of a geodesic \eqn{\gamma(t)} for \eqn{t \in (0,1)} which 
#' assumes two endpoints \eqn{\gamma(0)=X_1} and \eqn{\gamma(1)=X_2}.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{2} manifold-valued data where the first object is the starting point.
#' @param t a scalar in \eqn{(0,1)} for which the interpolation is taken.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return an interpolated object in matrix representation on \eqn{\mathcal{M}}.
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
#' ## FIND THE INTERPOLATED POINT AT "t=0.25"
#' mid.int = as.vector(riem.interp(sp.data, t=0.25, geometry="intrinsic"))
#' mid.ext = as.vector(riem.interp(sp.data, t=0.25, geometry="extrinsic"))
#' 
#' ## VISUALIZE
#' #  Prepare Lines and Points
#' thetas  = seq(from=0, to=pi/2, length.out=100)
#' quarter = cbind(cos(thetas), sin(thetas))
#' pic.pts = rbind(sp.start, mid.int, mid.ext, sp.end)
#' pic.col = c("black","red","green","black")
#' 
#' # Draw
#' opar <- par(no.readonly=TRUE)
#' par(pty="s")
#' plot(quarter, main="two interpolated points at t=0.25",
#'      xlab="x", ylab="y", type="l")
#' points(pic.pts, col=pic.col, pch=19)
#' text(mid.int[1]-0.1, mid.int[2], "intrinsic", col="red")
#' text(mid.ext[1]-0.1, mid.ext[2], "extrinsic", col="green")
#' par(opar)
#' 
#' @concept basic
#' @export
riem.interp <- function(riemobj, t=0.5, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.interp : input ",DNAME," should be an object of 'riemdata' class."))
  }
  if (length(riemobj$data)!=2){
    stop("* riem.interp : for 'riem.interp', please have two objects only.")
  }
  if ((length(t)>1)||(any(t<=0))||(any(t>=1))){
    stop("* riem.interp : 't' should be a single number in (0,1).")
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTATION
  vect   = rep(t, 2)
  out2   = basic_interpolate(riemobj$name, mygeom, riemobj$data[[1]], riemobj$data[[2]], vect)
  outmat = as.matrix(out2[,,1])
  return(outmat)
}