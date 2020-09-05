#' Fréchet Median and Variation
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' compute Fréchet median and variation with respect to the geometry by minimizing
#' \deqn{\textrm{min}_x \sum_{n=1}^N w_n \rho (x, x_n),\quad x\in\mathcal{M}} where
#' \eqn{\rho (x, y)} is a distance for two points \eqn{x,y\in\mathcal{M}}. 
#' If non-uniform weights are given, normalized version of the mean is computed 
#' and if \code{weight=NULL}, it automatically sets equal weights for all observations.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param weight weight of observations; if \code{NULL} it assumes equal weights, or a nonnegative length-\eqn{N} vector that sums to 1 should be given.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param maxiter maximum number of iterations to be run.
#' @param eps tolerance level for stopping criterion.
#' 
#' @return a named list containing\describe{
#' \item{median}{a median matrix on \eqn{\mathcal{M}}.}
#' \item{variation}{sum of (weighted) distances.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #        Example on Sphere : points near (0,1) on S^1 in R^2
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' ndata = 50
#' mydat = array(0,c(ndata,2))
#' for (i in 1:ndata){
#'   tgt = c(stats::rnorm(1, sd=2), 1)
#'   mydat[i,] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydat)
#' 
#' ## COMPUTE TWO MEANS
#' med.int = as.vector(riem.median(myriem, geometry="intrinsic")$median)
#' med.ext = as.vector(riem.median(myriem, geometry="extrinsic")$median)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(mydat[,1], mydat[,2], pch=19, xlim=c(-1.1,1.1), ylim=c(0,1.1),
#'      main="BLUE-extrinsic vs RED-intrinsic")
#' arrows(x0=0,y0=0,x1=med.int[1],y1=med.int[2],col="red")
#' arrows(x0=0,y0=0,x1=med.ext[1],y1=med.ext[2],col="blue")
#' par(opar)
#' 
#' @concept inference
#' @export
riem.median <- function(riemobj, weight=NULL, geometry=c("intrinsic","extrinsic"), 
                        maxiter=50, eps=1e-5){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.median : input ",DNAME," should be an object of 'riemdata' class."))
  }
  N = length(riemobj$data)
  if ((length(weight)==0)&&(is.null(weight))){
    myweight = rep(1/N, N)
  } else {
    myweight = check_weight(weight, N, "riem.median")
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  myiter = max(50, round(maxiter))
  myeps  = min(max(as.double(eps),0),1e-5)
  
  ## MAIN COMPUTATION
  if (all(mygeom=="intrinsic")){
    output = inference_median_intrinsic(riemobj$name, riemobj$data, myweight, myiter, myeps)
  } else {
    output = inference_median_extrinsic(riemobj$name, riemobj$data, myweight, myiter, myeps)
  }
  
  ## WRAP AND RETURN
  output$distvec = NULL # remove distance vector
  return(output)
}
