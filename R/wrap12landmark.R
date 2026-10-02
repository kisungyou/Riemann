#' Wrap Landmark Data on Shape Space
#' 
#' One of the frameworks used in shape space is to represent the data as landmarks. 
#' Each shape is a point set of \eqn{k} points in \eqn{\mathbf{R}^p} where each 
#' point is a labeled object. We consider general landmarks in \eqn{p=2,3,\ldots}. 
#' Note that when \eqn{p > 2}, it is stratified space but we assume singularities do not exist or 
#' are omitted. The wrapper takes translation and scaling out from the data to make it 
#' \emph{preshape} (centered, unit-norm). Also, for convenience, orthogonal 
#' Procrustes analysis is applied with the first observation being the reference so 
#' that all the other data are aligned to that reference. The full orthogonal
#' group O(p), including reflections, is used. This is not an SO(p)-only
#' orientation-preserving analysis. Singular configurations are rejected.
#' Fitted tangent models retain and apply their own training reference.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(k\times p\times n)} array where each slice along 3rd dimension is a \eqn{k}-ad in \eqn{\mathbf{R}^p}.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{k}-ads.}
#' }
#' 
#' @param reference optional reference configuration for alignment. If NULL,
#'   the first observation is used.
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of preshapes in \eqn{\mathbf{R}^p}.}
#'   \item{size}{size of each preshape.}
#'   \item{name}{name of the manifold of interests, \emph{"landmark"}}
#' }
#' 
#' @examples
#' ## USE 'GORILLA' DATA
#' data(gorilla)
#' riemobj = wrap.landmark(gorilla$male)
#' 
#' @references 
#' \insertRef{dryden_statistical_2016}{Riemann}
#' 
#' @concept wrapper
#' @export
wrap.landmark <- function(input, reference = NULL) {
  if (is.array(input)) {
    if (!check_3darray(input, symmcheck = FALSE)) {
      stop("Landmark input must be a nonempty three-dimensional array.", call. = FALSE)
    }
    data <- lapply(seq_len(dim(input)[3L]), function(i) {
      matrix(input[, , i, drop = FALSE], nrow = dim(input)[1L],
             dimnames = dimnames(input)[1:2])
    })
  } else if (is.list(input)) data <- input
  else stop("Landmark input must be a list of matrices or a three-dimensional array.", call. = FALSE)
  if (!check_list_eqsize(data) || !all(vapply(data, is.matrix, logical(1)))) {
    stop("Landmark configurations must be nonempty matrices of the same size.", call. = FALSE)
  }
  if (ncol(data[[1L]]) < 2L || nrow(data[[1L]]) <= ncol(data[[1L]])) {
    stop("Regular landmark configurations require at least p+1 landmarks in p >= 2 dimensions.", call. = FALSE)
  }
  nm <- dimnames(data[[1L]])
  if (!all(vapply(data, function(x) identical(dimnames(x), nm), logical(1)))) {
    stop("Landmark and coordinate names must have the same order in every observation.", call. = FALSE)
  }
  data <- lapply(data, aux_landmark_nearest)
  reference <- if (is.null(reference)) data[[1L]] else aux_landmark_nearest(reference)
  if (!identical(dim(reference), dim(data[[1L]]))) stop("Reference has incompatible dimensions.", call. = FALSE)
  data <- lapply(data, function(x) aux_landmark_match(reference, x))
  data <- lapply(data, function(x) { dimnames(x) <- nm; x })
  structure(list(data = data, size = dim(data[[1L]]), name = "landmark",
    dimnames = nm, representation = "centered_unit_preshape_O(p)",
    alignment_reference = reference, preprocessing = list(center = TRUE, scale = TRUE, reflections = TRUE)),
    class = "riemdata")
}

#' @keywords internal
#' @noRd
aux_landmark_nearest <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || !length(x) || any(!is.finite(x))) {
    stop("A landmark configuration must be a finite real matrix.", call. = FALSE)
  }
  scale <- max(abs(x))
  if (scale == 0) stop("A constant landmark configuration has no shape.", call. = FALSE)
  y <- sweep(x / scale, 2L, colMeans(x / scale), "-")
  size <- sqrt(sum(y * y))
  if (size == 0) stop("A constant landmark configuration has no shape.", call. = FALSE)
  y <- y / size
  values <- svd(y, nu = 0L, nv = 0L)$d
  if (length(values) < ncol(y) || min(values) <= 64 * .Machine$double.eps * max(values)) {
    stop("Singular landmark configurations are outside the supported regular shape space.", call. = FALSE)
  }
  y
}

#' @keywords internal
#' @noRd
aux_landmark_match <- function(x, y) {
  decomposition <- svd(crossprod(x, y))
  y %*% decomposition$v %*% t(decomposition$u)
}
