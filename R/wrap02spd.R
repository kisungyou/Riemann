#' Prepare Data on Symmetric Positive-Definite (SPD) Manifold
#' 
#' The collection of symmetric positive-definite matrices is a well-known example 
#' of matrix manifold. It is defined as
#' \deqn{\mathcal{S}_{++}^p = \lbrace X \in \mathbf{R}^{p\times p} ~\vert~ X^\top = X,~ \textrm{rank}(X)=p \rbrace}
#' where the rank condition means it is strictly positive definite. Please note that 
#' the geometry involving semi-definite matrices is considered in \code{wrap.spdk}. 
#' 
#' @param input SPD data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{an \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a SPD matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} SPD matrices.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times p)} correlation matrices.}
#'   \item{size}{size of each correlation matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"spd"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  Generate 5 observations; empirical covariance of normal observations.
#' #-------------------------------------------------------------------
#' #  Data Generation
#' d1 = array(0,c(3,3,5))
#' d2 = list()
#' for (i in 1:5){
#'   dat = matrix(rnorm(10*3),ncol=3)
#'   d1[,,i] = stats::cov(dat)
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.spd(d1)
#' test2 = wrap.spd(d2)
#' 
#' @concept wrapper
#' @export
wrap.spd <- function(input){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=TRUE)){
      stop("* wrap.spd : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.spd : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=TRUE)){
    stop("* wrap.spd : elements are not of same size.")
  }
  #  3. check
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = check_spd(tmpdata[[n]], n)
  }  
  
  # WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "spd"
  return(structure(output, class="riemdata"))
}
#' @keywords internal
#' @noRd
check_spd <- function(x, id){
  p = nrow(x)
  cond1 = (nrow(x)==ncol(x))
  cond2 = (round(mat_rank(x))==p) # full-rank
  cond3 = isSymmetric(x)
  if (cond1&&cond2&&cond3){
    return(x)
  } else {
    remainder = (id%%10)
    if (remainder==1){
      stop(paste0(" wrap.spd : ",id,"st object is not a valid SPD object."))
    } else if (remainder==2){
      stop(paste0(" wrap.spd : ",id,"nd object is not a valid SPD object."))
    } else if (remainder==3){
      stop(paste0(" wrap.spd : ",id,"rd object is not a valid SPD object."))
    } else {
      stop(paste0(" wrap.spd : ",id,"th object is not a valid SPD object."))
    }
  }
}