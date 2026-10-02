#' Prepare Data on Multinomial Manifold
#' 
#' Multinomial manifold is referred to the strictly positive data that sums to 1.
#' Also known as probability simplex or positive orthant, we denote \eqn{(p-1)} simplex 
#' in \eqn{\mathbf{R}^p} by 
#' \deqn{\Delta^{p-1} = \lbrace
#' x \in \mathbf{R}^p~\vert~ \sum_{i=1}^p x_i = 1, x_i > 0
#' \rbrace}
#' in that data are positive \eqn{L_1} unit-norm vectors. 
#' In \code{wrap.multinomial}, normalization is applied when each data point is not on the simplex, 
#' accepting positive counts or masses. Zero, negative, and nonfinite inputs are
#' rejected before normalization; no pseudocount is added.
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
wrap.multinomial <- function(input) {
  data <- riem_vector_input(input, "multinomial")
  data <- lapply(seq_along(data), function(i) {
    x <- data[[i]]
    x[] <- single_multinomial(as.numeric(x), i)
    x
  })
  riem_wrap_matrices(data, "multinomial")
}

single_multinomial <- function(vec, id) {
  if (!is.numeric(vec) || is.complex(vec) || length(vec) < 2L ||
      any(!is.finite(vec)) || any(vec <= 0)) {
    stop("Observation ", id, " must have at least two finite strictly positive entries.", call. = FALSE)
  }
  output <- vec / max(vec)
  output <- output / sum(output)
  if (any(output <= 0) || any(output >= 1)) {
    stop("Observation ", id, " is too close to the simplex boundary for finite precision.", call. = FALSE)
  }
  output
}
