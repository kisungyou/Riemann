#' Prepare Data on (Compact) Stiefel Manifold
#' 
#' Stiefel manifold \eqn{St(k,p)} is the set of \eqn{k}-frames in \eqn{\mathbf{R}^p}, 
#' which is indeed a Riemannian manifold. For usage in \pkg{Riemann} package, 
#' each data point is represented as a matrix by the convention
#' \deqn{St(k,p) = \lbrace X \in \mathbf{R}^{p\times k} ~\vert~ X^\top X = I_k \rbrace}
#' which means that columns are orthonormal. When the provided matrix is not 
#' an orthonormal basis as above, \code{wrap.stiefel} applies orthogonalization 
#' to extract valid basis information.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(p\times k\times n)} array where each slice along 3rd dimension is a \eqn{k}-frame.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times k)} \eqn{k}-frames.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{k}-frame orthonormal matrices.}
#'   \item{size}{size of each \eqn{k}-frame basis matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"stiefel"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  Generate 5 observations in St(2,4)
#' #-------------------------------------------------------------------
#' #  Data Generation by QR Decomposition
#' d1 = array(0,c(4,2,5))
#' d2 = list()
#' for (i in 1:5){
#'   d1[,,i] = qr.Q(qr(matrix(rnorm(4*2),ncol=2)))
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.stiefel(d1)
#' test2 = wrap.stiefel(d2)
#' 
#' @concept wrapper
#' @export
wrap.stiefel <- function(input){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=FALSE)){
      stop("* wrap.stiefel : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.stiefel : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=FALSE)){
    stop("* wrap.stiefel : elements are not of same size.")
  }
  #  3. check and transform to Stiefel
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = check_stiefel(tmpdata[[n]])
  }
  
  ## WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "stiefel"
  return(structure(output, class="riemdata"))
}
#' @keywords internal
#' @noRd
check_stiefel <- function(mat){
  k   = ncol(mat)
  tgt = t(mat)%*%mat
  eps = (base::norm(tgt-diag(k),"F")/sqrt(k))
  if (eps > 1e-10){
    return(base::qr.Q(base::qr(mat)))
  } else {
    return(mat)
  }
}
