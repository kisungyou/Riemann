#' Riemannian Manifold Metric Learning
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}} and 
#' corresponding label information, \code{riem.rmml} computes pairwise distance of data under Riemannian Manifold Metric Learning 
#' (RMML) framework based on equivariant embedding. When the number of data points 
#' is not sufficient, an inverse of scatter matrix does not exist analytically so 
#' the small regularization parameter \eqn{\lambda} is recommended with default value of \eqn{\lambda=0.1}.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param label a length-\eqn{N} vector of class labels. \code{NA} values are omitted.
#' @param lambda regularization parameter. If \eqn{\lambda \leq 0}, no regularization is applied.
#' @param as.dist logical; if \code{TRUE}, it returns \code{dist} object, else it returns a symmetric matrix.
#' 
#' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #            Distance between Two Classes of SPD Matrices
#' #
#' #  Class 1 : Empirical Covariance from Standard Normal Distribution
#' #  Class 2 : Empirical Covariance from Perturbed 'iris' dataset
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' data(iris)
#' ndata  = 10
#' mydata = list()
#' for (i in 1:ndata){
#'   mydata[[i]] = stats::cov(matrix(rnorm(100*4),ncol=4))
#' }
#' for (i in (ndata+1):(2*ndata)){
#'   tmpdata = as.matrix(iris[,1:4]) + matrix(rnorm(150*4,sd=0.5),ncol=4)
#'   mydata[[i]] = stats::cov(tmpdata)
#' }
#' myriem = wrap.spd(mydata)
#' mylabs = rep(c(1,2), each=ndata)
#' 
#' ## COMPUTE GEODESIC AND RMML PAIRWISE DISTANCE
#' pdgeo = riem.pdist(myriem)
#' pdmdl = riem.rmml(myriem, label=mylabs)
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(pdgeo[,(2*ndata):1], main="geodesic distance", axes=FALSE)
#' image(pdmdl[,(2*ndata):1], main="RMML distance", axes=FALSE)
#' par(opar)
#' 
#' @references 
#' \insertRef{zhu_generalized_2018}{Riemann}
#' 
#' @concept learning
#' @export
riem.rmml <- function(riemobj, label, lambda=0.1, as.dist=FALSE){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.rmml : input ",DNAME," should be an object of 'riemdata' class."))
  }
  mylbd   = max(0, as.double(lambda))
  mydist  = as.logical(as.dist)
  mylabel = as.factor(label)
  
  if (any(is.na(mylabel))){
    idselect = (!is.na(mylabel))
    mylabel  = mylabel[idselect]
    mydata   = riemobj$data[idselect]
  } else {
    mydata   = riemobj$data
  }
  mylabel = as.integer(mylabel)-1
  
  ## COMPUTE
  distmat = learning_rmml(riemobj$name, mydata, mylbd, mylabel)
  if (mydist){
    return(stats::as.dist(distmat))
  } else {
    return(distmat)
  }
}