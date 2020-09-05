#' Prepare Data on Multinomial Manifold
#' 
#' Multinomial manifold is referred to the data that is nonnegative and sums to 1. 
#' Also known as probability simplex or positive orthant, we denote \eqn{(p-1)} simplex 
#' in \eqn{\mathbf{R}^p} by 
#' \deqn{\Delta^{p-1} = \lbrace
#' x \in \mathbf{R}^p~\vert~ \sum_{i=1}^p x_i = 1, x_i > 0
#' \rbrace}
#' in that data are positive \eqn{L_1} unit-norm vectors. 
#' In \code{wrap.multinomial}, normalization is applied when each data point is not on the simplex, 
#' but if vectors contain values not in \eqn{(0,1)}, it returns errors.
#' 
#' @param input data vectors to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{matrix}{an \eqn{(n \times p)} matrix of row observations.}
#' \item{list}{a length-\eqn{n} list whose elements are length-\eqn{p} vectors.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times 1)} matrices in \eqn{\Delta^{p-1}}.}
#'   \item{size}{dimension of the ambient space.}
#'   \item{name}{name of the manifold of interests, \emph{"multinomial"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' d1 = array(0,c(5,3))
#' d2 = list()
#' for (i in 1:5){
#'   single  = abs(stats::rnorm(3))
#'   d1[i,]  = single
#'   d2[[i]] = single
#' }
#' 
#' ## RUN
#' test1 = wrap.multinomial(d1)
#' test2 = wrap.multinomial(d2)
#' 
#' @concept wrapper
#' @export
wrap.multinomial <- function(input){
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
    stop("* wrap.multinomial : input should be either a 2d matrix or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=FALSE)){
    stop("* wrap.multinomial : elements are not vectors of same size.")
  }
  #  3. check each element
  N = length(tmpdata)
  for (n in 1:N){
    tgtvec = single_multinomial(as.vector(tmpdata[[n]]), n)
    tmpdata[[n]] = matrix(tgtvec, ncol = 1)
  }
  
  ## WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "multinomial"
  return(structure(output, class="riemdata"))
}
#' @keywords internal
#' @noRd
single_multinomial <- function(vec, id){
  output = vec/base::sum(vec)
  if (any(output <= 0)||any(output >= 1)){
    stop(paste0("* wrap.multinomial : ",id,"-th vector is not a suitable object. See the description."))
  }
  return(output)
}