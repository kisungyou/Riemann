#' Kernel Principal Component Analysis
#' 
#' Although the method of Kernel Principal Component Analysis (KPCA) was originally 
#' developed to visualize non-linearly distributed data in Euclidean space, 
#' we graft this to the case for manifolds where extrinsic geometry is explicitly available.
#' The algorithm uses Gaussian kernel with 
#' \deqn{K(X_i, X_j) = \exp\left( - \frac{d^2 (X_i, X_j)}{2 \sigma^2} \right )}
#' where \eqn{\sigma} is a bandwidth parameter and \eqn{d(\cdot, \cdot)} is 
#' an extrinsic distance defined on a specific manifold. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param sigma the bandwidth parameter (default: 1).
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' \item{vars}{a length-\eqn{N} vector of eigenvalues from kernelized covariance matrix.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example for Gorilla Skull Data : 'gorilla'
#' #-------------------------------------------------------------------
#' ## PREPARE THE DATA
#' #  Aggregate two classes into one set
#' data(gorilla)
#' 
#' mygorilla = array(0,c(8,2,59))
#' for (i in 1:29){
#'   mygorilla[,,i] = gorilla$male[,,i]
#' }
#' for (i in 30:59){
#'   mygorilla[,,i] = gorilla$female[,,i-29]
#' }
#' 
#' gor.riem = wrap.landmark(mygorilla)
#' gor.labs = c(rep("red",29), rep("blue",30))
#' 
#' ## APPLY KPCA WITH DIFFERENT KERNEL BANDWIDTHS
#' kpca1 = riem.kpca(gor.riem, sigma=0.01)
#' kpca2 = riem.kpca(gor.riem, sigma=1)
#' kpca3 = riem.kpca(gor.riem, sigma=100)

#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(kpca1$embed, pch=19, col=gor.labs, main="sigma=1/100")
#' plot(kpca2$embed, pch=19, col=gor.labs, main="sigma=1")
#' plot(kpca3$embed, pch=19, col=gor.labs, main="sigma=100")
#' par(opar)
#' 
#' @references 
#' \insertRef{scholkopf_kernel_1997}{Riemann}
#' 
#' @concept visualization
#' @export
riem.kpca <- function(riemobj, ndim=2, sigma=1.0){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.kpca : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim  = max(2, round(ndim))
  mysigma = max(sqrt(.Machine$double.eps), as.double(sigma))
  
  ## COMPUTE VIA RCPP AND RETURN
  output = visualize_kpca(riemobj$name, riemobj$data, mysigma, myndim)
  output$vars = as.vector(output$vars)
  return(output)
}