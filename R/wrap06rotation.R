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
wrap.rotation <- function(input) {
  data <- riem_matrix_input(input, square = TRUE)
  data <- lapply(data, function(x) { single_rotcheck(x); x })
  riem_wrap_matrices(data, "rotation")
}

single_rotcheck <- function(x, id = 0L) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || any(!is.finite(x)) ||
      !nrow(x) || nrow(x) != ncol(x) ||
      norm(crossprod(x) - diag(nrow(x)), "F") > 1e-10 * sqrt(nrow(x)) ||
      abs(det(x) - 1) > 1e-10) {
    stop("Observation ", id, " must be an orthogonal rotation matrix with determinant +1.", call. = FALSE)
  }
  invisible(TRUE)
}
