#' Prepare Data on Correlation Manifold
#' 
#' The collection of correlation matrices is considered as a subset (and quotient) of 
#' the well-known SPD manifold. In our package, it is defined as
#' \deqn{\mathcal{C}_{++}^p = \lbrace X \in \mathbf{R}^{p\times p} ~\vert~ X^\top = X,~ X \succ 0,~ \textrm{diag}(X) = 1 \rbrace}
#' where \eqn{X \succ 0} means strictly positive definite. Please note that
#' the geometry involving semi-definite correlation matrices is not the objective here. 
#' 
#' @param input correlation data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{an \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a correlation matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} correlation matrices.}
#' }
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times p)} correlation matrices.}
#'   \item{size}{size of each correlation matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"correlation"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #
#' #  5 observations; empirical correlation of normal observations.
#' #-------------------------------------------------------------------
#' #  Data Generation
#' d1 = array(0,c(3,3,5))
#' d2 = list()
#' for (i in 1:5){
#'   dat = matrix(rnorm(10*3),ncol=3)
#'   d1[,,i] = stats::cor(dat)
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.correlation(d1)
#' test2 = wrap.correlation(d2)
#' 
#' @concept wrapper
#' @export
wrap.correlation <- function(input) {
  data <- riem_matrix_input(input, square = TRUE)
  data <- lapply(data, function(x) check_corr(x, 0L))
  riem_wrap_matrices(data, "correlation")
}

check_corr <- function(x, id) {
  if (!check_spdmat(x) || any(abs(diag(x) - 1) > 64 * .Machine$double.eps)) {
    stop("Observation ", id, " must be symmetric positive definite with unit diagonal.", call. = FALSE)
  }
  x / 2 + t(x) / 2
}
