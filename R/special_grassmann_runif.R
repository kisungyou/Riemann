#' Generate Uniform Samples on Grassmann Manifold
#' 
#' It generates \eqn{n} random samples from Grassmann manifold \eqn{Gr(k,p)}.
#' 
#' @param n number of samples to be generated.
#' @param k dimension of the subspace.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of \eqn{(p\times k)} basis of \eqn{k}-subspaces.}
#' \item{\code{"array"}}{a \eqn{(p\times k\times n)} 3D array whose slices are \eqn{k}-subspace basis.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.grassmann}} for more details.}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' @seealso \code{\link{stiefel.runif}}, \code{\link{wrap.grassmann}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Draw Samples on Grassmann Manifold 
#' #-------------------------------------------------------------------
#' #  Multiple Return Types with 3 Observations of 5-dim subspaces in R^10
#' dat.list = grassmann.runif(n=3, k=5, p=10, type="list")
#' dat.arr3 = grassmann.runif(n=3, k=5, p=10, type="array")
#' dat.riem = grassmann.runif(n=3, k=5, p=10, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @concept grassmann
#' @export
grassmann.runif <- function(n, k, p, type=c("list","array","riemdata")){
  #  PREPROCESSING
  n = round(n)
  p = round(p) # k-frame in R^p
  k = round(k)
  if (k > p){
    stop("* grassmann.runif : 'k <= p' is a required condition.")
  }
  retype = ifelse(missing(type),"riemdata",
                  match.arg(tolower(type), c("list","array","riemdata")))
  
  #  GENERATE, WRAP, RETURN
  tmpout = runif_stiefel(p,k,n) # C++ version
  if (all(retype=="array")){
    return(tmpout)
  } else {
    tmpobj = wrap.grassmann(tmpout)
    if (all(retype=="list")){
      return(tmpobj$data)
    } else {
      return(tmpobj)
    }
  }
}