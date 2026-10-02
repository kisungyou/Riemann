#' Prepare Data on Symmetric Positive-Definite (SPD) Manifold
#' 
#' The collection of symmetric positive-definite matrices is a well-known example 
#' of matrix manifold. It is defined as
#' \deqn{\mathcal{S}_{++}^p = \lbrace X \in \mathbf{R}^{p\times p} ~\vert~ X^\top = X,~ v^\top Xv>0~\textrm{for all}~v\ne0 \rbrace}
#' where positivity is verified by Cholesky factorization after a scale-relative symmetry check.
#' No ridge or positive-definite projection is applied. Tiny physical units alone
#' do not make a matrix ill-conditioned. Please note that
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
#'   \item{data}{a list of \eqn{(p\times p)} SPD matrices.}
#'   \item{size}{size of each SPD matrix.}
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
wrap.spd <- function(input) {
  if (is.array(input)) {
    if (!check_3darray(input, symmcheck = TRUE)) {
      stop("wrap.spd requires a nonempty square three-dimensional array.", call. = FALSE)
    }
    tmpdata <- lapply(seq_len(dim(input)[3L]), function(i) {
      matrix(input[, , i, drop = FALSE], nrow = dim(input)[1L],
             dimnames = dimnames(input)[1:2])
    })
  } else if (is.list(input)) {
    tmpdata <- input
  } else stop("wrap.spd requires a list of matrices or a three-dimensional array.", call. = FALSE)
  if (!check_list_eqsize(tmpdata, check.square = TRUE)) {
    stop("wrap.spd requires nonempty square matrices of the same dimensions.", call. = FALSE)
  }
  tmpdata <- lapply(seq_along(tmpdata), function(i) check_spd(tmpdata[[i]], i))
  if (!all(vapply(tmpdata, function(x) identical(dimnames(x), dimnames(tmpdata[[1L]])), logical(1)))) {
    stop("wrap.spd requires consistent feature ordering, including whether names are supplied.", call. = FALSE)
  }
  structure(list(data = tmpdata, size = dim(tmpdata[[1L]]), name = "spd",
                 dimnames = dimnames(tmpdata[[1L]])), class = "riemdata")
}

#' @keywords internal
#' @noRd
check_spd <- function(x, id) {
  if (!check_spdmat(x)) {
    stop("wrap.spd: observation ", id,
         " must be a finite real symmetric positive-definite matrix.", call. = FALSE)
  }
  # Only remove symmetry error already bounded by the validation tolerance.
  x / 2 + t(x) / 2
}
