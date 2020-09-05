#' Prepare Data on Grassmann Manifold
#' 
#' Grassmann manifold \eqn{Gr(k,p)} is the set of \eqn{k}-planes, or \eqn{k}-dimensional subspaces in \eqn{R^p}, 
#' which means that for a given matrix \eqn{Y \in \mathbf{R}{p\times k}}, the column space \eqn{SPAN(Y)} is an element 
#' in Grassmann manifold. We use a convention that each element in \eqn{Gr(k,p)} is represented as an orthonormal basis (ONB) \eqn{X \in \mathbf{R}^{p\times k}} where
#' \deqn{X^\top X = I_k.} If not provided in such a form, this wrapper takes a QR decomposition of the given data 
#' to recover a corresponding ONB.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{an \eqn{(p\times k\times n)} array where each slice along 3rd dimension is a \eqn{k}-subspace basis in dimension \eqn{p}.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times k)} basis for \eqn{k}-subspace.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of k-subspace basis matrices.}
#'   \item{size}{size of each k-subspace basis matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"grassmann"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  Generate 5 observations in Gr(2,4)
#' #-------------------------------------------------------------------
#' #  Generation
#' d1 = array(0,c(4,2,5))
#' d2 = list()
#' for (i in 1:5){
#'   d1[,,i] = matrix(rnorm(4*2), ncol=2)
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.grassmann(d1)
#' test2 = wrap.grassmann(d2)
#' 
#' @concept wrapper
#' @export
wrap.grassmann <- function(input){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=FALSE)){
      stop("* wrap.grassmann : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.grassmann : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=FALSE)){
    stop("* wrap.grassmann : elements are not of same size.")
  }
  #  3. check and transform to Stiefel
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = check_stiefel(tmpdata[[n]])
  }
  
  ############################################################
  # WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "grassmann"
  return(structure(output, class="riemdata"))
}