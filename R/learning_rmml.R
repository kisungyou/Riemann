#' Riemannian Manifold Metric Learning
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}} and 
#' corresponding label information, \code{riem.rmml} computes pairwise distance of data under Riemannian Manifold Metric Learning 
#' (RMML) framework based on equivariant embedding. When the number of data points 
#' is not sufficient, an inverse of scatter matrix does not exist analytically so 
#' the small regularization parameter \eqn{\lambda} is recommended with default value of \eqn{\lambda=0.1}.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param label a length-\eqn{N} vector of class labels. \code{NA} values are omitted.
#' @param lambda regularization parameter. If \eqn{\lambda \leq 0}, no regularization is applied.
#' @param as.dist logical; if \code{TRUE}, it returns \code{dist} object, else it returns a symmetric matrix.
#' 
#' @return a S3 \code{dist} object or symmetric matrix of pairwise distances for
#' the observations with nonmissing labels, in their original order, according
#' to \code{as.dist}. At least two labeled observations are required.
#' @details Embedding dimensions are determined from the manifold's actual
#' equivariant embedding. In particular, a Grassmann frame is represented by its
#' full orthogonal projector, so the result is unchanged by a change of basis
#' within each represented subspace.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #            Distance between Two Classes of SPD Matrices
#' #
#' #  Class 1 : Empirical Covariance from Standard Normal Distribution
#' #  Class 2 : Empirical Covariance from Perturbed 'iris' dataset
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' data(iris)
#' ndata  = 10
#' mydata = list()
#' for (i in 1:ndata){
#'   mydata[[i]] = stats::cov(matrix(rnorm(100*4),ncol=4))
#' }
#' for (i in (ndata+1):(2*ndata)){
#'   tmpdata = as.matrix(iris[,1:4]) + matrix(rnorm(150*4,sd=0.5),ncol=4)
#'   # This matrix-entry demonstration uses positional coordinates in both groups.
#'   mydata[[i]] = unname(stats::cov(tmpdata))
#' }
#' myriem = wrap.spd(mydata)
#' mylabs = rep(c(1,2), each=ndata)
#' 
#' ## COMPUTE GEODESIC AND RMML PAIRWISE DISTANCE
#' pdgeo = riem.pdist(myriem)
#' pdmdl = riem.rmml(myriem, label=mylabs)
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(pdgeo[,(2*ndata):1], main="geodesic distance", axes=FALSE)
#' image(pdmdl[,(2*ndata):1], main="RMML distance", axes=FALSE)
#' par(opar)
#' 
#' @references 
#' \insertRef{zhu_generalized_2018}{Riemann}
#' 
#' @concept learning
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
riem.rmml <- function(riemobj, label, lambda=0.1, as.dist=FALSE) {
  riem_validate_data(riemobj)
  if (!is.atomic(label) || !is.null(dim(label)) || length(label) != length(riemobj$data)) {
    stop("label must be a vector with one label per observation.", call. = FALSE)
  }
  if (!is.numeric(lambda) || is.complex(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("lambda must be a finite real number.", call. = FALSE)
  }
  if (!is.logical(as.dist) || length(as.dist) != 1L || is.na(as.dist)) {
    stop("as.dist must be TRUE or FALSE.", call. = FALSE)
  }
  keep <- !is.na(label)
  if (sum(keep) < 2L) stop("At least two nonmissing labels are required.", call. = FALSE)
  labels <- as.integer(factor(label[keep])) - 1L
  distances <- learning_rmml(riemobj$name, riemobj$data[keep], max(0, lambda), labels)
  if (as.dist) stats::as.dist(distances) else distances
}
