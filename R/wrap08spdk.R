#' Prepare Data on Positive Semidefinite Manifold of Fixed Rank
#' 
#' When \eqn{(p\times p)} positive semidefinite matrices are of fixed rank \eqn{k < p}, they form
#' a geometric structure represented by \eqn{(p\times k)} matrices,
#' \deqn{SPD(k,p) = \lbrace X \in \mathbf{R}^{(p\times p)}~\vert~ Y Y^\top = X, \textrm{rank}(X) = k \rbrace}
#' It's key difference from \eqn{\mathcal{S}_{++}^p} is that all matrices should be 
#' of fixed rank \eqn{k} where \eqn{k} is usually smaller than \eqn{p}. Inputs are 
#' given as \eqn{(p\times p)} matrices with specified \eqn{k} and \code{wrap.spdk} 
#' decomposes matrices whose numerical rank is exactly \eqn{k}. Higher-rank inputs
#' are rejected; rank reduction must be performed explicitly before wrapping.
#' Numerical rank uses eigenvalues exceeding \eqn{64p} times machine precision
#' relative to the largest absolute eigenvalue. Negative eigenvalues beyond this
#' roundoff tolerance are rejected.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a rank-\eqn{k} matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} matrices of rank-\eqn{k}.}
#' }
#' @param k rank of the positive semidefinite matrices, an integer in \eqn{[1,p-1]}.
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times k)} representation of the corresponding rank-\eqn{k} SPSD matrices.}
#'   \item{size}{size of each representation matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"spdk"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #-------------------------------------------------------------------
#' #  Data Generation
#' d1 = array(0,c(10,10,3))
#' d2 = list()
#' for (i in 1:3){
#'   dat = matrix(rnorm(10*2),nrow=10)
#'   d1[,,i] = tcrossprod(dat)
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.spdk(d1, k=2)
#' test2 = wrap.spdk(d2, k=2)
#' 
#' 
#' @references 
#' \insertRef{journee_lowrank_2010}{Riemann}
#' 
#' @concept wrapper
#' @export
wrap.spdk <- function(input, k) {
  data <- riem_matrix_input(input, square = TRUE)
  if (!is.numeric(k) || is.complex(k) || length(k) != 1L || !is.finite(k) ||
      k != floor(k) || k < 1L || k >= nrow(data[[1L]])) {
    stop("k must be one integer between 1 and p-1 for the fixed-rank semidefinite geometry.", call. = FALSE)
  }
  data <- lapply(seq_along(data), function(i) single_spdkcheck(data[[i]], i, as.integer(k)))
  result <- riem_wrap_matrices(data, "spdk")
  result$representation <- "full_rank_factor_modulo_O(k)"
  result
}

single_spdkcheck <- function(x, id, k) {
  scale <- max(abs(x))
  tol <- 64 * .Machine$double.eps * nrow(x)
  if (scale == 0 || max(abs(x / scale - t(x / scale))) > tol) {
    stop("Observation ", id, " must be a nonzero symmetric positive semidefinite matrix.", call. = FALSE)
  }
  eig <- eigen((x / scale) / 2 + t(x / scale) / 2, symmetric = TRUE)
  cutoff <- tol * max(abs(eig$values))
  if (min(eig$values) < -cutoff || sum(eig$values > cutoff) != k) {
    stop("Observation ", id, " must be positive semidefinite with numerical rank exactly k; no truncation is applied.", call. = FALSE)
  }
  result <- sweep(eig$vectors[, seq_len(k), drop = FALSE], 2L,
                  sqrt(eig$values[seq_len(k)]) * sqrt(scale), "*")
  rownames(result) <- rownames(x)
  result
}
