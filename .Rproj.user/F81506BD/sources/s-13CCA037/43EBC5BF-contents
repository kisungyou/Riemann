#' Principal Geodesic Analysis
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' Principal Geodesic Analysis (PGA) finds a low-dimensional embedding by decomposing 
#' 2nd-order information in tangent space at an intrinsic mean of the data. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim an integer-valued target dimension.
#' 
#' @return a named list containing \describe{
#' \item{center}{an intrinsic mean in a matrix representation form.}
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # 10 perturbed data points near (0,0,1) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata = list()
#' for (i in 1:10){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 21:30){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' mylabs = rep(c(1,2,3), each=10)
#' 
#' ## EMBEDDING WITH MDS AND PGA
#' embed2mds = riem.mds(myriem, ndim=2, geometry="intrinsic")$embed
#' embed2pga = riem.pga(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2mds, main="Multidimensional Scaling",    col=mylabs, pch=19)
#' plot(embed2pga, main="Principal Geodesic Analysis", col=mylabs, pch=19)
#' par(opar)
#' 
#' @references 
#' \insertRef{fletcher_principal_2004}{Riemann}
#' 
#' @concept visualization
#' @export
riem.pga <- function(riemobj, ndim=2){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.pga : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myndim = max(2, round(ndim))
  
  ## PRELIMINARY COMPUTE
  tmprun = visualize_pga(riemobj$name, riemobj$data)
  return(tmprun)
  # nrows  = dim(tmprun$logvec)[1]
  # ncols  = dim(tmprun$logvec)[2]
  # N      = dim(tmprun$logvec)[3]
  # 
  # ## PCA
  # rowlogs = array(0,c(N,nrows*ncols))
  # for (n in 1:N){
  #   rowlogs[n,] = as.vector(tmprun$logvec[,,n])
  # }
  # embed   = rowlogs%*%base::eigen(t(rowlogs)%*%rowlogs)$vectors[,1:myndim]
  # 
  # ## WRAP AND RETURN
  # tmprun$logvec = NULL
  # tmprun$embed  = embed
  # return(tmprun)
}