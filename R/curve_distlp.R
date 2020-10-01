#' Distance between Two Curves on Manifolds
#' 
#' Given two curves \eqn{\gamma_1, \gamma_2 : I \rightarrow \mathcal{M}}, we are 
#' interested in measuring the discrepancy of two curves. Usually, data are given 
#' as discrete observations so we are offering several methods to perform the task. See 
#' the section below for detailed description.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data along the curve.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data along the curve.
#' @param vect a vector of domain values. If given \code{Null} (default), sequence \code{1:N} is set.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{p}{an exponent (default: 2).}
#' }
#' 
#' @return the distance value.
#' 
#' @section Default Method : Trapezoidal Approximation
#' Assume \eqn{\gamma_1 (t_i) = X_i} and \eqn{\gamma_2 (t_i) = Y_i} for 
#' \eqn{i=1,2,\ldots,N}. In the Euclidean space, \eqn{L_p} distance between two 
#' scalar-valued functions is defined as 
#' \deqn{L_p^p (\gamma_1 (x), \gamma_2 (x) = \int_{\mathcal{X}} |\gamma_1 (x) - \gamma_2 (x)|^p dx }. 
#' We extend this approach to manifold-valued curves
#' \deqn{L_p^p (\gamma_1 (t), \gamma_2 (t)) = \int_{t\in I} d^p (\gamma_1 (t), \gamma_2 (t)) dt}
#' where \eqn{d(\cdot,\cdot)} is an intrinsic/extrinsic distance on manifolds. With the given 
#' representations, the above integral is approximated using trapezoidal rule.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                          Curves on Sphere
#' #
#' #  curve1 : y = 0.5*cos(x) on the tangent space at (0,0,1)
#' #  curve2 : y = 0.5*cos(x) on the tangent space at (0,0,1)
#' #  curve3 : y = 0.5*sin(x) on the tangent space at (0,0,1)
#' #
#' # * distance between curve1 & curve2 should be close to 0.
#' # * distance between curve1 & curve3 should be large.
#' #-------------------------------------------------------------------
#' ## GENERATION
#' vecx  = seq(from=-0.9, to=0.9, length.out=50)
#' vecy1 = 0.5*cos(vecx) + rnorm(50, sd=0.05)
#' vecy2 = 0.5*cos(vecx) + rnorm(50, sd=0.05)
#' vecy3 = 0.5*sin(vecx) + rnorm(50, sd=0.05)
#' 
#' ## WRAP AS RIEMOBJ
#' mat1 = cbind(vecx, vecy1, 1); mat1 = mat1/sqrt(rowSums(mat1^2))
#' mat2 = cbind(vecx, vecy2, 1); mat2 = mat2/sqrt(rowSums(mat2^2))
#' mat3 = cbind(vecx, vecy3, 1); mat3 = mat3/sqrt(rowSums(mat3^2))
#' 
#' rcurve1 = wrap.sphere(mat1)
#' rcurve2 = wrap.sphere(mat2)
#' rcurve3 = wrap.sphere(mat3)
#' 
#' ## COMPUTE DISTANCES
#' riem.distlp(rcurve1, rcurve2, vect=vecx)
#' riem.distlp(rcurve1, rcurve3, vect=vecx)
#' 
#' @concept curve
#' @export
riem.distlp <- function(riemobj1, riemobj2, vect=NULL, geometry=c("intrinsic","extrinsic"), ...){
  ## PREPARE : EXPLICIT
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* riem.distlp : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* riem.distlp : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  if (!all(riemobj1$name==riemobj2$name)){
    stop("* riem.distlp : two inputs are from different manifolds.")
  }
  if (!all(riemobj1$size==riemobj2$size)){
    stop("* riem.distlp : two inputs are of different size.")
  }
  if (length(riemobj1$data)!=length(riemobj2$data)){
    stop("* riem.distlp : two inputs should be of same size.")
  }
  N = length(riemobj1$data)
  if ((length(vect)==0)&&(is.null(vect))){
    myvect = 1:N
  } else {
    myvect = as.vector(vect)
  }
  if (length(myvect)!=N){
    stop("* riem.distlp : length of 'vect' is not matching to that of the data.")
  }
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## PREPARE : IMPLICIT
  pars     = list(...)
  pnames   = names(pars)
  myp      = max(1, ifelse(("p"%in%pnames), as.double(pars$p), 2))
  mymethod = "lp"
  
  # COMPUTE THE DISTANCE
  outdist = switch(mymethod,
                   lp = curvedist_lp(riemobj1$name, mygeometry, riemobj1$data, riemobj2$data, myvect, myp))
  return(outdist)
}
