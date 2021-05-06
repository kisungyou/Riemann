#' Manifold-to-Scalar Kernel Regression with K-Fold Cross Validation 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data corresponding to \eqn{X_1,\ldots,X_N}.
#' @param y a length-\eqn{N} vector of dependent variable values.
#' @param bandwidths a vector of nonnegative numbers that control smoothness.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param kfold the number of folds for cross validation.
#' 
#' @return a named list of S3 class \code{m2skreg} containing
#' \describe{
#' \item{ypred}{a length-\eqn{N} vector of optimal smoothed responses.}
#' \item{bandwidth}{the optimal bandwidth value.}
#' \item{inputs}{a list containing both \code{riemobj} and \code{y} for future use.}
#' \item{errors}{a matrix whose columns are \code{bandwidths} values and corresponding errors measure in SSE.}
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
#' set.seed(496) 
#' npts = 100
#' nlev = 0.25
#' thetas = seq(from=0, to=pi/2, length.out=npts)
#' Xstack = cbind(rep(0,npts), sin(thetas), cos(thetas))
#' 
#' Xriem  = wrap.sphere(Xstack)
#' ytrue  = sin(seq(from=0, to=2*pi, length.out=npts))
#' ynoise = ytrue + rnorm(npts, sd=nlev)
#' 
#' # FIT WITH 5-FOLD CV
#' cv_band = (10^seq(from=-4, to=-1, length.out=200))
#' cv_fit  = riem.m2skregCV(Xriem, ynoise, bandwidths=cv_band)
#' cv_err  = cv_fit$errors
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(1:npts, cv_fit$ypred, pch=19, cex=0.5, "b", xlab="", main="optimal prediction")
#' lines(1:npts, ytrue, col="red", lwd=1.5)
#' plot(cv_err[,1], cv_err[,2], "b", pch=19, cex=0.5, main="5-fold CV errors",
#'      xlab="bandwidth", ylab="SSE")
#' abline(v=cv_fit$bandwidth, col="blue", lwd=1.5)
#' par(opar)
#' }
#' 
#' @concept inference
#' @export
riem.m2skregCV <- function(riemobj, y, bandwidths=seq(from=0.01, to=1, length.out=10), geometry=c("intrinsic","extrinsic"), kfold=5){
  # CHECK INPUTS
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.m2skregCV : input ",DNAME," should be an object of 'riemdata' class."))
  }
  N = length(riemobj$data)
  y = as.vector(y)
  if (length(y)!=N){
    stop(paste0("* riem.m2skregCV : length of 'y' should equal to ",N,"."))
  }
  mybandvecs = base::pmax(sqrt(.Machine$double.eps), as.double(bandwidths))
  mygeometry = ifelse(missing(geometry),"intrinsic",
                       match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  mynfolds  = max(2, round(kfold))
  mydistmat = basic_pdist(riemobj$name, riemobj$data, mygeometry)

  # CV
  splitgp  <- get_splits(N, mynfolds)
  myerrors <- rep(0, length(mybandvecs))
  for (i in 1:length(mybandvecs)){
    # current bandwidth parameter
    now_bandwidth = mybandvecs[i]
    # variable for saving the errors
    now_error = 0
    for (j in 1:length(mynfolds)){
      # separate data
      id_j  = splitgp[[j]]
      # use the auxiliary function to compute SSE
      # sse_j = riem.m2skregCV.single(riemobj$name, riemobj$data, y, id_j, now_bandwidth, mygeometry)
      sse_j = riem.m2skregCV.each(mydistmat, y, id_j, now_bandwidth)
      # update the error term
      now_error = now_error + sse_j
    }
    myerrors[i] = now_error
  }
  
  # What is the Optimal Bandwidth
  opt.bandwidth = as.vector(mybandvecs[which.min(myerrors)])[1]
  opt.ypred     = riem.m2skreg(riemobj, y, bandwidth=opt.bandwidth, geometry=mygeometry)$ypred
  
  # record the errors & remove any
  errormat = cbind(mybandvecs, myerrors)
  colnames(errormat) = c("bandwidth","SSE")
  
  idremove = which(is.na(myerrors))
  errormat = errormat[-idremove,]
  
  ## RETURN THE OUTPUT
  output = list()
  output$ypred = opt.ypred
  output$bandwidth = opt.bandwidth
  output$inputs = list(riemobj, y)
  output$errors = errormat
  return(structure(output, class="m2skreg"))
}


#' @keywords internal
#' @noRd
get_splits <- function(N, K){
  return(suppressWarnings(split(sample(1:N, N), as.factor(1:K))))
}
#' @keywords internal
#' @noRd
riem.m2skregCV.each <- function(pdistmat, y, id.now, bandwidth){
  # parameters
  M = length(id.now)
  N = length(y)-M
  
  # separate Y's
  train_y = as.vector(y[-id.now])
  test_y  = as.vector(y[id.now])
  
  # compute
  distmat = pdistmat[-id.now, id.now]
  pred_y  = rep(0,M)
  for (m in 1:M){
    tgtvec = as.vector(distmat[,m])
    tgtscd = base::exp(-(tgtvec^2)/(2*(bandwidth^2)))
    pred_y[m] = base::sum(tgtscd*train_y)/base::sum(tgtscd)
  }
  
  # Return SSE
  return(base::sum((pred_y - test_y)^2))
}
#' #' @keywords internal
#' #' @noRd
#' riem.m2skregCV.single <- function(X.name, X.data, y.data, id.now, bandwidth, geometry){
#'   # separate the data
#'   train_X = X.data[-id.now]
#'   train_y = as.vector(y.data[-id.now])
#'   
#'   test_X = X.data[id.now]
#'   test_y = y.data[id.now]
#'   
#'   # information
#'   N = length(train_y)
#'   M = length(test_y)
#'   
#'   # Pairwise Distance matrix of size (N x M) 
#'   distmat = basic_pdist2(X.name, train_X, test_X, geometry)
#'   
#'   # Do the prediction
#'   pred_y = rep(0,M)
#'   for (m in 1:M){
#'     tgtvec = as.vector(distmat[,m])
#'     tgtscd = base::exp(-(tgtvec^2)/(2*(bandwidth^2)))
#'     pred_y[m] = base::sum(tgtscd*train_y)/base::sum(tgtscd)
#'   }
#'   
#'   # Return SSE
#'   return(base::sum((pred_y - test_y)^2))
#' }


