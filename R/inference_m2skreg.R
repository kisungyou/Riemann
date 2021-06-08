#' Manifold-to-Scalar Kernel Regression 
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}} and 
#' scalars \eqn{y_1, y_2, \ldots, y_N \in \mathbf{R}}, perform the Nadaraya-Watson kernel 
#' regression by 
#' \deqn{\hat{m}_h (X) = \frac{\sum_{i=1}^n K \left( \frac{d(X,X_i)}{h}  \right) y_i}{\sum_{i=1}^n K \left( \frac{d(X,X_i)}{h}  \right)}}
#' where the Gaussian kernel is defined as
#' \deqn{K(x) := \frac{1}{\sqrt{2\pi}} \exp \left( - \frac{x^2}{2}\right)} 
#' with the bandwidth parameter \eqn{h > 0} that controls the degree of smoothness. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data corresponding to \eqn{X_1,\ldots,X_N}.
#' @param y a length-\eqn{N} vector of dependent variable values.
#' @param bandwidth a nonnegative number that controls smoothness.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' 
#' @return a named list of S3 class \code{m2skreg} containing
#' \describe{
#' \item{ypred}{a length-\eqn{N} vector of smoothed responses.}
#' \item{bandwidth}{the bandwidth value that was originally provided, which is saved for future use.}
#' \item{inputs}{a list containing both \code{riemobj} and \code{y} for future use.}
#' }
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #                    Example on Sphere S^2
#' #
#' #  X : equi-spaced points from (0,0,1) to (0,1,0)
#' #  y : sin(x) with perturbation
#' #-------------------------------------------------------------------
#' # GENERATE DATA
#' npts = 100
#' nlev = 0.25
#' thetas = seq(from=0, to=pi/2, length.out=npts)
#' Xstack = cbind(rep(0,npts), sin(thetas), cos(thetas))
#' 
#' Xriem  = wrap.sphere(Xstack)
#' ytrue  = sin(seq(from=0, to=2*pi, length.out=npts))
#' ynoise = ytrue + rnorm(npts, sd=nlev)
#' 
#' # FIT WITH DIFFERENT BANDWIDTHS
#' fit1 = riem.m2skreg(Xriem, ynoise, bandwidth=0.001)
#' fit2 = riem.m2skreg(Xriem, ynoise, bandwidth=0.01)
#' fit3 = riem.m2skreg(Xriem, ynoise, bandwidth=0.1)
#' 
#' # VISUALIZE
#' xgrd <- 1:npts
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(xgrd, fit1$ypred, pch=19, cex=0.5, "b", xlab="", ylim=c(-2,2), main="h=1e-3")
#' lines(xgrd, ytrue, col="red", lwd=1.5)
#' plot(xgrd, fit2$ypred, pch=19, cex=0.5, "b", xlab="", ylim=c(-2,2), main="h=1e-2")
#' lines(xgrd, ytrue, col="red", lwd=1.5)
#' plot(xgrd, fit3$ypred, pch=19, cex=0.5, "b", xlab="", ylim=c(-2,2), main="h=1e-1")
#' lines(xgrd, ytrue, col="red", lwd=1.5)
#' par(opar)
#' }
#' 
#' @concept inference
#' @export
riem.m2skreg <- function(riemobj, y, bandwidth=0.5, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.m2skreg : input ",DNAME," should be an object of 'riemdata' class."))
  }
  N = length(riemobj$data)
  y = as.vector(y)
  if (length(y)!=N){
    stop(paste0("* riem.m2skreg : length of 'y' should equal to ",N,"."))
  }
  mybandwidth = max(sqrt(.Machine$double.eps), as.double(bandwidth))
  mygeometry  = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## COMPUTE PAIRWISE DISTANCE
  distmat  = basic_pdist(riemobj$name, riemobj$data, mygeometry)
  
  ## FILL IN THE DIFFERENCES
  ypred = rep(0,N)
  for (n in 1:N){
    tgtvec = as.vector(distmat[n,])
    tgtscd = exp(-(tgtvec^2)/(2*(mybandwidth^2)))
    ypred[n] = base::sum(tgtscd*y)/base::sum(tgtscd)
  }
  
  ## RETURN THE OUTPUT
  output = list()
  output$ypred = ypred
  output$bandwidth = mybandwidth
  output$inputs = list(riemobj, y)
  return(structure(output, class="m2skreg"))
}



#' Prediction for Manifold-to-Scalar Kernel Regression 
#' 
#' Given new observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}}, plug in 
#' the data with respect to the fitted model for prediction. 
#' 
#' @param object an object of \code{m2skreg} class. See \code{\link{riem.m2skreg}} for more details.
#' @param newdata a S3 \code{"riemdata"} class for manifold-valued data corresponding to \eqn{X_1,\ldots,X_M}.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return a length-\eqn{M} vector of predictted values.
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #                    Example on Sphere S^2
#' #
#' #  X : equi-spaced points from (0,0,1) to (0,1,0)
#' #  y : sin(x) with perturbation
#' #
#' #  Our goal is to check whether the predict function works well
#' #  by comparing the originally predicted values vs. those of the same data.
#' #-------------------------------------------------------------------
#' # GENERATE DATA
#' npts = 100
#' nlev = 0.25
#' thetas = seq(from=0, to=pi/2, length.out=npts)
#' Xstack = cbind(rep(0,npts), sin(thetas), cos(thetas))
#' 
#' Xriem  = wrap.sphere(Xstack)
#' ytrue  = sin(seq(from=0, to=2*pi, length.out=npts))
#' ynoise = ytrue + rnorm(npts, sd=nlev)
#' 
#' # FIT & PREDICT
#' obj_fit   = riem.m2skreg(Xriem, ynoise, bandwidth=0.01)
#' yval_fits = obj_fit$ypred
#' yval_pred = predict(obj_fit, Xriem)
#' 
#' # VISUALIZE
#' xgrd <- 1:npts
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(xgrd, yval_fits, pch=19, cex=0.5, "b", xlab="", ylim=c(-2,2), main="original fit")
#' lines(xgrd, ytrue, col="red", lwd=1.5)
#' plot(xgrd, yval_pred, pch=19, cex=0.5, "b", xlab="", ylim=c(-2,2), main="from 'predict'")
#' lines(xgrd, ytrue, col="red", lwd=1.5)
#' par(opar)
#' }
#' 
#' @seealso \code{\link{riem.m2skreg}}
#' @concept inference
#' @export
predict.m2skreg <- function(object, newdata, geometry=c("intrinsic","extrinsic"), ...){
  # Check Inputs
  if (!inherits(object,"m2skreg")){
    stop("* predict : input is not an object of 'm2skreg' class.")
  }
  mygeometry  = ifelse(missing(geometry),"intrinsic",
                       match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  my_old  = object$inputs[[1]] # original X's
  my_yvec = object$inputs[[2]] # original y's
  my_bdh  = as.double(object$bandwidth)
  
  my_new  = newdata
  if (!check_tworiems(my_old, my_new)){
    stop("* predict : input 'newdata' is a not valid object due to one of many possible reasons.")
  }
  N = length(my_old$data) # reference points
  M = length(my_new$data) # provided data length
  
  # Pairwise Distance matrix of size (N x M) 
  distmat = basic_pdist2(my_old$name, my_old$data, my_new$data, mygeometry)
  
  # Do the prediction
  ypred = rep(0,M)
  for (m in 1:M){
    tgtvec = as.vector(distmat[,m])
    tgtscd = base::exp(-(tgtvec^2)/(2*(my_bdh^2)))
    ypred[m] = base::sum(tgtscd*my_yvec)/base::sum(tgtscd)
  }
  
  # Return
  return(ypred)
}