#' Build Lightweight Coreset
#' 
#' Given manifold-valued data \eqn{X_1,X_2,\ldots,X_N \in \mathcal{M}}, this algorithm 
#' draws a weighted coreset using the lightweight importance-sampling scheme
#' proposed by the reference below. The Euclidean approximation theorem is not
#' asserted for arbitrary manifold geometries.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param M positive integer number of independent draws (default: \eqn{\lceil N/2 \rceil}). Values greater than \eqn{N} are permitted.
#' @param geometry (case-insensitive) name or saved specification of a geometry supporting means. Legacy aliases \code{"intrinsic"} and \code{"extrinsic"} are accepted; see \code{\link{riem.geometry}}.
#' @param ... extra parameters including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{coreid}{a length-\eqn{M} index vector; repeated indices are retained.}
#' \item{weight}{a length-\eqn{M} vector of importance weights, \eqn{1/(M q_i)} for each drawn index \eqn{i}.}
#' }
#' 
#' @details Each index is drawn independently with replacement, with probability
#' \eqn{q_i = 1/(2N) + d(X_i,\mu)^2/(2\sum_j d(X_j,\mu)^2)}, where
#' \eqn{\mu} is the computed mean in the selected geometry. When all distances
#' are zero the probabilities are uniform. Thus, for any fixed centers, the
#' importance-weighted coreset cost is an unbiased estimate of the full-data
#' squared-distance cost. This identity does not require the computed mean to
#' be globally optimal. Duplicate draws count separately and must not be
#' discarded without summing their weights. Use \code{set.seed()} for reproducibility.
#'
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # * 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # * 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # * 10 perturbed data points near (0,0,1) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata = list()
#' for (i in 1:10){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 21:30){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' 
#' ## MDS FOR VISUALIZATION
#' embed2 = riem.mds(myriem, ndim=2)$embed
#' 
#' ## FIND CORESET OF SIZES 3, 6, 9
#' core1 = riem.coreset18B(myriem, M=3)
#' core2 = riem.coreset18B(myriem, M=6)
#' core3 = riem.coreset18B(myriem, M=9)
#' 
#' col1 = rep(1,30); col1[core1$coreid] = 2
#' col2 = rep(1,30); col2[core2$coreid] = 2
#' col3 = rep(1,30); col3[core3$coreid] = 2
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embed2, pch=19, col=col1, main="coreset size=3")
#' plot(embed2, pch=19, col=col2, main="coreset size=6")
#' plot(embed2, pch=19, col=col3, main="coreset size=9")
#' par(opar)
#' 
#' @references 
#' \insertRef{bachem_scalable_2018a}{Riemann}
#' 
#' @concept learning
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
riem.coreset18B <- function(riemobj, M=max(1L, ceiling(length(riemobj$data)/2)),
                           geometry=c("intrinsic","extrinsic"), ...) {
  riem_validate_data(riemobj)
  if (missing(geometry)) geometry <- NULL
  geometry <- riem_resolve_geometry(riemobj, geometry, capability = "mean")
  M <- riem_kmeans_integer(M, "M")
  pars <- riem_legacy_parameters(list(...), c("maxiter", "eps"))
  maxiter <- if (is.null(pars$maxiter)) 50L else riem_kmeans_integer(pars$maxiter, "maxiter")
  eps <- if (is.null(pars$eps)) 1e-5 else riem_legacy_positive(pars$eps, "eps")
  result <- learning_coreset18B(riemobj$name, geometry$backend,
                              riemobj$data, M, maxiter, eps)
  indices <- as.integer(result$id) + 1L
  list(coreid = indices, weight = 1 / (M * as.vector(result$qx)[indices]))
}
