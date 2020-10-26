#' Generate Random Samples from Multivariate Normal Distribution
#' 
#' In \eqn{\mathbf{R}^p}, random samples are drawn
#' \deqn{X_1,X_2,\ldots,X_n~ \sim ~ \mathcal{N}(\mu, \Sigma)}
#' where \eqn{\mu \in \mathbf{R}^p} is a mean vector and \eqn{\Sigma \in \textrm{SPD}(p)} 
#' is a positive definite covariance matrix.
#' 
#' @param n the number of samples to be generated.
#' @param mu mean vector.
#' @param sigma covariance matrix.
#' 
#' @return either (1) a length-\eqn{p} vector (\eqn{n=1}) or (2) an \eqn{(n\times p)} matrix where rows are random samples.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #   Generate Random Data and Compare with Empirical Covariances
#' #
#' # In R^5 with zero mean and diagonal covariance, 
#' # generate 100 and 200 observations and compute MLE covariance.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mymu  = rep(0,5)
#' mysig = diag(5)
#' 
#' ## MLE FOR COVARIANCE
#' smat1 = stats::cov(rmvnorm(n=100, mymu, mysig))
#' smat2 = stats::cov(rmvnorm(n=200, mymu, mysig))
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(mysig[,5:1], axes=FALSE, main="true covariance")
#' image(smat1[,5:1], axes=FALSE, main="empirical cov with n=100")
#' image(smat2[,5:1], axes=FALSE, main="empirical cov with n=200")
#' par(opar)
#' 
#' @concept utility
#' @export
rmvnorm <- function(n=1, mu, sigma){
  myn   = max(round(n), 1)
  mymu  = as.vector(mu)
  mysig = as.matrix(sigma)
  
  output = cpp_rmvnorm(myn, mymu, mysig)
  if (myn < 2){
    return(as.vector(output))
  } else {
    return(output)
  }
}