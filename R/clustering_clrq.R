#' Competitive Learning Riemannian Quantization
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform clustering via Competitive Learning Riemannian Quantization (CLRQ). 
#' Originally, the algorithm is designed for finding voronoi cells that are 
#' used in domain quantization. Given the discrete measure of data, centers of the cells 
#' play a role of cluster centers and data are labeled accordingly based on the distance 
#' to voronoi centers. For an iterative update of centers, gradient descent algorithm 
#' adapted for the Riemannian manifold setting is used with the gain factor sequence
#' \deqn{\gamma_t = \frac{a}{1 + b \sqrt{t}}}
#' where two parameters \eqn{a,b} are represented by \code{par.a} and \code{par.b}. For 
#' initialization, we provide k-means++ and random seeding options as in k-means.
#' 
#' @seealso \code{\link{riem.kmeans}}
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param k the number of clusters.
#' @param init (case-insensitive) name of an initialization scheme. (default: \code{"plus"}.)
#' @param gain.a parameter \eqn{a} for gain factor sequence.
#' @param gain.b parameter \eqn{b} for gain factor sequence.
#' 
#' @return a named list containing\describe{
#' \item{centers}{a 3d array where each slice along 3rd dimension is a matrix representation of class centers.}
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # class 1 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # class 3 : 10 perturbed data points near (0,0,1) on S^2 in R^3
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
#' ## CLRQ WITH K=2,3,4
#' clust2 = riem.clrq(myriem, k=2)
#' clust3 = riem.clrq(myriem, k=3)
#' clust4 = riem.clrq(myriem, k=4)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(mds2d, pch=19, main="true label", col=mylabs)
#' plot(mds2d, pch=19, main="K=2", col=clust2$cluster)
#' plot(mds2d, pch=19, main="K=3", col=clust3$cluster)
#' plot(mds2d, pch=19, main="K=4", col=clust4$cluster)
#' par(opar)
#' 
#' @references 
#' \insertRef{lebrigant_quantization_2019}{Riemann}
#' 
#' \insertRef{bonnabel_stochastic_2013}{Riemann}
#' 
#' @concept clustering
#' @export
riem.clrq <- function(riemobj, k=2, init=c("plus","random"), gain.a = 1, gain.b = 1){
  ## PREPARE
  N          = length(riemobj$data)
  par.init   = ifelse(missing(init), "plus", match.arg(tolower(init),c("plus","random")))
  par.k      = max(1, round(k))
  par.a      = max(sqrt(.Machine$double.eps), as.double(gain.a)) # gamma_t = a/(1 + b t^0.5) : Bonnabel's family; stochastic gradient descent
  par.b      = max(sqrt(.Machine$double.eps), as.double(gain.b))
  par.ab     = (par.a/(1+par.b)) 
  if ((par.ab <= 0)||(par.ab >= 1)){
    stop("* riem.clrq : for 'gain.a' and 'gain.b' parameters, 'gain.a/(1+gain.b)' should be a number in (0,1).")
  }
  
  ## INITIALIZATION
  if (all(par.init=="random")){
    rec.id0 = base::sample(1:N, k)
  } else {
    distobj     = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, "intrinsic"))
    func.import = utils::getFromNamespace("hidden_kmeanspp", "maotai")
    rec.id0     = as.vector(func.import(distobj, k=par.k)$center)
  }
  
  ## MAIN COMPUTATION
  output = clustering_clrq(riemobj$name, riemobj$data, rec.id0, par.a, par.b)
  
  ## WRAP AND RETURN
  output$cluster = base::apply(output$pdist2, 1, which.min)
  output$pdist2  = NULL
  return(output)
}