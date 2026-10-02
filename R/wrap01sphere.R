#' Prepare Data on Sphere
#' 
#' The unit hypersphere (sphere, for short) is one of the most fundamental curved 
#' space in studying geometry. Precisely, we denote \eqn{(p-1)} sphere in \eqn{\mathbf{R}^p} by
#' \deqn{\mathcal{S}^{p-1} = \lbrace x \in \mathbf{R}^p ~ \vert ~ x^\top x = \|x\|^2 = 1 \rbrace}
#' where vectors are of unit norm. In \code{wrap.sphere}, normalization is applied when 
#' each data point is not on the unit sphere.
#' 
#' @param input data vectors to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{matrix}{an \eqn{(n \times p)} matrix of row observations of unit norm.}
#' \item{list}{a length-\eqn{n} list whose elements are length-\eqn{p} vectors of unit norm.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times 1)} matrices in \eqn{\mathcal{S}^{p-1}}.}
#'   \item{size}{dimension of the ambient space.}
#'   \item{name}{name of the manifold of interests, \emph{"sphere"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  Generate 5 observations in S^2 embedded in R^3.
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
#' test1 = wrap.sphere(d1)
#' test2 = wrap.sphere(d2)
#' 
#' @concept wrapper
#' @export
wrap.sphere <- function(input) {
  tmpdata <- riem_vector_input(input, "sphere")
  structure(list(data = tmpdata, size = dim(tmpdata[[1L]]), name = "sphere",
                 dimnames = dimnames(tmpdata[[1L]])), class = "riemdata")
}
