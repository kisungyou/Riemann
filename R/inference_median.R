#' Fr\ifelse{html}{\out{&eacute;}}{\ifelse{latex}{\out{\'e}}{e}}chet Median and Variation
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' compute Fr\ifelse{html}{\out{&eacute;}}{\ifelse{latex}{\out{\'e}}{e}}chet median and variation with respect to the geometry by minimizing
#' \deqn{\textrm{min}_x \sum_{n=1}^N w_n \rho (x, x_n),\quad x\in\mathcal{M}} where
#' \eqn{\rho (x, y)} is a distance for two points \eqn{x,y\in\mathcal{M}}. 
#' If non-uniform weights are given, normalized version of the median is computed
#' and if \code{weight=NULL}, it automatically sets equal weights for all observations.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param weight Finite nonnegative observation weights with positive sum. They
#'   are normalized internally. \code{NULL} uses equal weights. Zero-weight
#'   observations are validated but do not enter the calculation.
#' @param geometry Geometry name or specification. \code{NULL} selects the
#'   default geometry. Legacy \code{"intrinsic"} and \code{"extrinsic"}
#'   aliases remain available for supported combinations. SPD names include
#'   \code{"affine_invariant"} and \code{"log_euclidean"}.
#' @param ... Named controls: \describe{
#'   \item{maxiter}{Positive maximum number of accepted iterations (default 50).}
#'   \item{eps}{Positive absolute tolerance for the compatible subgradient residual
#'     (default \code{1e-5}).}
#'   \item{max_backtrack}{Positive maximum number of trial steps per iteration,
#'     at most 1024 (default 50).}
#'   \item{trace}{Whether to retain the iteration history (default \code{FALSE}).}
#'   \item{init}{Optional initial matrix with the observation dimensions.
#'     Closed-form calculations do not use an initializer.}
#' } Unknown controls are errors.
#' 
#' @details Intrinsic medians use a modified Weiszfeld direction with
#'   backtracking. Coincident observations retain their subgradient mass, and
#'   an observation is accepted as a nonsmooth solution only after the relevant
#'   certificate is checked. On nonconvex manifolds, stationarity is local.
#'
#'   For the log-Euclidean SPD chart, the extrinsic route computes the geometric
#'   median of matrix logarithms and exponentiates it. For a curved embedding
#'   such as the sphere or Grassmann manifold, the returned estimator is the
#'   projection of an ambient geometric median. It need not minimize the sum
#'   of chordal distances constrained to the manifold. The \code{estimand} and
#'   \code{diagnostic_scope} fields identify this distinction; convergence and
#'   trace describe the ambient optimization. \code{ambient_objective} records
#'   its cost before projection. An ambiguous inverse projection is an error.
#'
#' @return A \code{riem_summary} object retaining \code{median} and
#'   \code{variation}, plus \code{objective}, resolved \code{geometry}, normalized
#'   \code{weights}, \code{converged}, \code{termination}, \code{iterations},
#'   \code{subgradient_residual} when applicable, controls, and an optional
#'   \code{trace}. \code{variation} equals the normalized weighted sum of
#'   distances at the returned matrix; for a projected ambient median this is
#'   a descriptive manifold cost, not the optimized ambient objective.
#'   An unsuccessful iteration warns and returns the last accepted estimate.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #        Example on Sphere : points near (0,1) on S^1 in R^2
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' ndata = 50
#' mydat = array(0,c(ndata,2))
#' for (i in 1:ndata){
#'   tgt = c(stats::rnorm(1, sd=2), 1)
#'   mydat[i,] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydat)
#' 
#' ## COMPUTE TWO MEANS
#' med.int = as.vector(riem.median(myriem, geometry="intrinsic")$median)
#' med.ext = as.vector(riem.median(myriem, geometry="extrinsic")$median)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(mydat[,1], mydat[,2], pch=19, xlim=c(-1.1,1.1), ylim=c(0,1.1),
#'      main="BLUE-extrinsic vs RED-intrinsic")
#' arrows(x0=0,y0=0,x1=med.int[1],y1=med.int[2],col="red")
#' arrows(x0=0,y0=0,x1=med.ext[1],y1=med.ext[2],col="blue")
#' par(opar)
#' 
#' @concept inference
#' @export
riem.median <- function(riemobj, weight = NULL, geometry = NULL, ...) {
  riem_fit_summary(riemobj, weight, geometry, list(...), "median", match.call())
}
