#' Prepare Data on Euclidean Space
#' 
#' Euclidean space \eqn{\mathbf{R}^p} is the most common space for data analysis, 
#' which can be considered as a Riemannian manifold with flat metric. Since
#' the space of matrices is isomorphic to Euclidean space after vectorization, 
#' we consider the inputs as \eqn{p}-dimensional vectors.
#' 
#' @param input data vectors to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{matrix}{an \eqn{(n \times p)} matrix of row observations.}
#' \item{list}{a length-\eqn{n} list whose elements are length-\eqn{p} vectors.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times 1)} matrices in \eqn{\mathbf{R}^p}.}
#'   \item{size}{dimension of the ambient space.}
#'   \item{name}{name of the manifold of interests, \emph{"euclidean"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  Generate 5 observations in R^3 in Matrix and List.
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' d1 = array(0,c(5,3))
#' d2 = list()
#' for (i in 1:5){
#'   single  = stats::rnorm(3)
#'   d1[i,]  = single
#'   d2[[i]] = single
#' }
#' 
#' ## RUN
#' test1 = wrap.euclidean(d1)
#' test2 = wrap.euclidean(d2)
#' 
#' @concept wrapper
#' @export
wrap.euclidean <- function(input){
  ## TAKE EITHER 2D ARRAY {n x p} OR A LIST
  #  1. data format
  if (is.matrix(input)){
    N = nrow(input)
    tmpdata = list()
    for (i in 1:N){
      tmpdata[[i]] = as.vector(input[i,])
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.euclidean : input should be either a 2d matrix or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=FALSE)){
    stop("* wrap.euclidean : elements are not vectors of same size.")
  }
  #  3. check each element
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = matrix(tmpdata[[n]], ncol=1)
  }
  
  ## WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "euclidean"
  return(structure(output, class="riemdata"))
}