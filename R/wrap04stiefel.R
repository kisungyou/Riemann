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
wrap.stiefel <- function(input) {
  data <- riem_matrix_input(input, square = FALSE)
  data <- lapply(data, function(x) check_stiefel(x))
  riem_wrap_matrices(data, "stiefel")
}

check_stiefel <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(mat) || is.complex(mat) ||
      any(!is.finite(mat)) || !riem_full_column_rank(mat)) {
    stop("A frame must be a finite real full-column-rank matrix with rows >= columns.", call. = FALSE)
  }
  if (norm(crossprod(mat) - diag(ncol(mat)), "F") <= 1e-10 * sqrt(ncol(mat))) return(mat)
  labels <- dimnames(mat)
  result <- qr.Q(qr(mat / max(abs(mat))))
  dimnames(result) <- labels
  result
}
