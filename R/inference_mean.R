#' Fr\ifelse{html}{\out{&eacute;}}{\ifelse{latex}{\out{\'e}}{e}}chet Mean and Variation
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' compute Fr\ifelse{html}{\out{&eacute;}}{\ifelse{latex}{\out{\'e}}{e}}chet mean and variation with respect to the geometry by minimizing
#' \deqn{\textrm{min}_x \sum_{n=1}^N w_n \rho^2 (x, x_n),\quad x\in\mathcal{M}} where
#' \eqn{\rho (x, y)} is a distance for two points \eqn{x,y\in\mathcal{M}}. 
#' If non-uniform weights are given, normalized version of the mean is computed 
#' and if \code{weight=NULL}, it automatically sets equal weights (\eqn{w_i = 1/n}) for all observations.
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
#'   \item{eps}{Positive absolute tolerance for the compatible gradient norm
#'     (default \code{1e-5}).}
#'   \item{max_backtrack}{Positive maximum number of trial steps per iteration,
#'     at most 1024 (default 50).}
#'   \item{trace}{Whether to retain the iteration history (default \code{FALSE}).}
#'   \item{init}{Optional initial matrix with the observation dimensions.
#'     Closed-form calculations do not use an initializer.}
#' } Unknown controls are errors.
#' 
#' @details Intrinsic iteration uses the average logarithm direction and an
#'   Armijo backtracking rule for the stated squared-distance objective.
#'   Euclidean and log-Euclidean SPD means use closed forms. On spaces with
#'   nonconvex objectives, termination with \code{"stationary"} establishes
#'   first-order stationarity to the requested tolerance only. A stationary
#'   point may be a local minimum, a saddle, or a local maximum; the criterion
#'   does not establish local minimality, uniqueness, or global optimality.
#'
#' @return A \code{riem_summary} object retaining \code{mean} and
#'   \code{variation}, plus \code{objective}, resolved \code{geometry}, normalized
#'   \code{weights}, \code{converged}, \code{termination}, \code{iterations},
#'   \code{gradient_norm}, \code{step_norm}, controls, and an optional \code{trace}.
#'   \code{variation} equals the normalized weighted sum of squared distances
#'   at the returned mean. An unsuccessful iteration warns and returns the last
#'   accepted estimate. For a curved extrinsic embedding, diagnostics and trace
#'   refer to the ambient calculation, as recorded in \code{diagnostic_scope}.
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
#' mean.int = as.vector(riem.mean(myriem, geometry="intrinsic")$mean)
#' mean.ext = as.vector(riem.mean(myriem, geometry="extrinsic")$mean)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(mydat[,1], mydat[,2], pch=19, xlim=c(-1.1,1.1), ylim=c(0,1.1),
#'      main="BLUE-extrinsic vs RED-intrinsic")
#' arrows(x0=0,y0=0,x1=mean.int[1],y1=mean.int[2],col="red")
#' arrows(x0=0,y0=0,x1=mean.ext[1],y1=mean.ext[2],col="blue")
#' par(opar)
#' 
#' @concept inference
#' @export
riem.mean <- function(riemobj, weight = NULL, geometry = NULL, ...) {
  riem_fit_summary(riemobj, weight, geometry, list(...), "mean", match.call())
}
