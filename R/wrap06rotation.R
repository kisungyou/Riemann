#' Prepare Data on Rotation Group
#' 
#' Rotation group, also known as special orthogonal group, is a Riemannian 
#' manifold
#' \deqn{SO(p) = \lbrace Q \in \mathbf{R}^{p\times p}~\vert~ Q^\top Q = I, \textrm{det}(Q)=1 \rbrace }
#' where the name originates from an observation that when \eqn{p=2,3} these matrices are rotation of 
#' shapes/configurations. 
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a rotation matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} rotation matrices.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times p)} rotation matrices.}
#'   \item{size}{size of each rotation matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"rotation"}}
#' }
#' 
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' d1 = array(0,c(3,3,5))
#' d2 = list()
#' for (i in 1:5){
#'   single  = qr.Q(qr(matrix(rnorm(9),nrow=3)))
#'   d1[,,i] = single
#'   d2[[i]] = single
#' }
#' 
#' ## RUN
#' test1 = wrap.rotation(d1)
#' test2 = wrap.rotation(d2)
#' 
#' @concept wrapper
#' @export
wrap.rotation <- function(input){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=FALSE)){
      stop("* wrap.rotation : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.rotation : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=TRUE)){
    stop("* wrap.rotation : elements are not of same size.")
  }
  #  3. check and transform to Stiefel
  N = length(tmpdata)
  for (n in 1:N){
    tmpcheck     = single_rotcheck(tmpdata[[n]], n)
    tmpdata[[n]] = tmpdata[[n]]
  }
  
  ## WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "rotation"
  return(structure(output, class="riemdata")) 
}
#' @keywords internal
#' @noRd
single_rotcheck <- function(x, id=0){
  p = nrow(x)
  if (nrow(x)!=ncol(x)){
    stop(paste0("* wrap.rotation : ",id,"-th element is not a square matrix."))
  }
  if ((norm((t(x)%*%x)-diag(p),"F")/sqrt(p) >= 1e-10)){
    stop(paste0("* wrap.rotation : ",id,"-th element does not satisfy X'*X = I."))
  }
  if (abs(base::det(x)-1) >= 1e-10){
    stop(paste0("* wrap.rotation : ",id,"-th element's determinant is not close to 1."))
  }
  return(TRUE)
}