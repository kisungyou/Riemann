#' Wrap Landmark Data on Shape Space
#' 
#' One of the frameworks used in shape space is to represent the data as landmarks. 
#' Each shape is a point set of \eqn{k} points in \eqn{\mathbf{R}^p} where each 
#' point is a labeled object. We consider general landmarks in \eqn{p=2,3,\ldots}. 
#' Note that when \eqn{p > 2}, it is stratified space but we assume singularities do not exist or 
#' are omitted. The wrapper takes translation and scaling out from the data to make it 
#' \emph{preshape} (centered, unit-norm). Also, for convenience, orthogonal 
#' Procrustes analysis is applied with the first observation being the reference so 
#' that all the other data are rotated to match the shape of the first.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(k\times p\times n)} array where each slice along 3rd dimension is a \eqn{k}-ad in \eqn{\mathbf{R}^p}.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{k}-ads.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of preshapes in \eqn{\mathbf{R}^p}.}
#'   \item{size}{size of each preshape.}
#'   \item{name}{name of the manifold of interests, \emph{"landmark"}}
#' }
#' 
#' @examples
#' ## USE 'GORILLA' DATA
#' data(gorilla)
#' riemobj = wrap.landmark(gorilla$male)
#' 
#' @references 
#' \insertRef{dryden_statistical_2016}{Riemann}
#' 
#' @concept wrapper
#' @export
wrap.landmark <- function(input){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=FALSE)){
      stop("* wrap.landmark : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.landmark : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=FALSE)){
    stop("* wrap.landmark : elements are not of same size.")
  }
  #  3. normalize
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = aux_landmark_nearest(tmpdata[[n]])
  }
  #  4. compute extrinsic mean
  rot.center = tmpdata[[1]]%*%base::solve(base::eigen(stats::cov(tmpdata[[1]]))$vectors)
  # rot.center = tmpdata[[1]]
  for (n in 1:N){
    tmpdata[[n]] = aux_landmark_match(rot.center, tmpdata[[n]])
  }
  
  ############################################################
  # WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "landmark"
  return(structure(output, class="riemdata"))
}


# auxiliary function for the wrapper --------------------------------------
#' @keywords internal
#' @noRd
aux_landmark_nearest <- function(x){
  y = x - matrix(rep(base::colMeans(x),nrow(x)),nrow=nrow(x),byrow=TRUE)
  return(y/base::norm(y,"F"))
}
#' @keywords internal
#' @noRd
aux_landmark_match <- function(x,y){
  sxy = base::svd(t(x)%*%y)
  return(y%*%sxy$v%*%t(sxy$u))
}