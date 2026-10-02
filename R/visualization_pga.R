#' Tangent Principal Component Analysis
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' Principal Geodesic Analysis (PGA) finds a low-dimensional embedding by decomposing 
#' 2nd-order information in tangent space at an intrinsic mean of the data. 
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param ndim positive integer requested number of components. Requests above
#'   numerical rank are truncated with a warning.
#' @param geometry a supported geometry name or saved geometry specification.
#'   Defaults to the geometry of the data. Supports affine-invariant and
#'   log-Euclidean SPD, round spheres, Euclidean observations, and regular
#'   orthogonal-quotient landmark shapes.
#' @param weight nonnegative observation weights, normalized to sum to one.
#' @param maxiter maximum iterations for the reference mean.
#' @param eps convergence tolerance for the reference mean.
#' @param rank.tol relative singular-value tolerance for numerical rank.
#' @param center.tangent whether to subtract the weighted tangent-coordinate
#'   mean before decomposition. The fitted offset is retained for reuse.
#'
#' @details This is tangent PCA, a local linear approximation, rather than an
#'   optimization over principal geodesic submanifolds. SPD coordinates preserve
#'   the selected metric using symmetric vectorization with square-root-of-two
#'   off-diagonal weights. Affine-invariant coordinates are whitened at the mean;
#'   log-Euclidean coordinates use the matrix-log chart. Sphere coordinates use
#'   an orthonormal tangent basis. Landmark coordinates use horizontal log
#'   vectors after alignment to the trained reference under the full orthogonal
#'   group (reflections included); singular alignments are rejected.
#'
#'   Component variances use normalized weights divided by
#'   \code{1 - sum(weight^2)}, agreeing with sample PCA for equal weights.
#'   This denominator is evaluated as twice the sum of positive pair products
#'   to avoid cancellation with highly unequal weights. Positive weights that
#'   underflow to zero during normalization are rejected. With only one
#'   positive weight, the centered fit has zero rank and zero variance by
#'   convention; this is not an estimate of a sample covariance. An uncentered
#'   fit requires at least two positive weights.
#'   Explained variance describes the tangent representation, not exact global
#'   manifold variation. A nonconverged reference mean causes an error.
#'   Prediction retains the trained mean, coordinate convention, offset and
#'   loadings. Reconstruction returns wrapped observations in this same frame.
#'   Sphere reconstruction is restricted to tangent norms below pi; landmark
#'   reconstruction is restricted below pi/2 and to regular shapes.
#' 
#' @return a named list containing \describe{
#' \item{center}{an intrinsic mean in a matrix representation form.}
#' \item{embed}{an \eqn{N}-by-effective-rank matrix of component scores.}
#' \item{loadings}{orthonormal directions in the stored tangent coordinates.}
#' \item{geometry}{the persistent geometry specification.}
#' \item{offset}{the fitted tangent-coordinate offset.}
#' \item{rank}{numerical rank of the weighted coordinate matrix.}
#' \item{variance}{component variances in decreasing order.}
#' \item{diagnostics}{reference-mean convergence and approximation information.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # 10 perturbed data points near (0,0,1) on S^2 in R^3
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
#' mylabs = rep(c(1,2,3), each=10)
#' 
#' ## EMBEDDING WITH MDS AND PGA
#' embed2mds = riem.mds(myriem, ndim=2, geometry="intrinsic")$embed
#' embed2pga = riem.pga(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2mds, main="Multidimensional Scaling",    col=mylabs, pch=19)
#' plot(embed2pga, main="Principal Geodesic Analysis", col=mylabs, pch=19)
#' par(opar)
#' 
#' @references 
#' \insertRef{fletcher_principal_2004}{Riemann}
#' 
#' @concept visualization
#' @export
riem.pga <- function(riemobj, ndim = 2, geometry = NULL, weight = NULL,
                     maxiter = 200, eps = 1e-8,
                     rank.tol = sqrt(.Machine$double.eps), center.tangent = TRUE) {
  riem_validate_data(riemobj)
  geom <- riem_resolve_geometry(riemobj, geometry, capability = "tangent")
  if (!is.numeric(ndim) || length(ndim) != 1L || !is.finite(ndim) ||
      ndim < 1 || ndim != floor(ndim)) {
    stop("ndim must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(rank.tol) || length(rank.tol) != 1L ||
      !is.finite(rank.tol) || rank.tol <= 0 || rank.tol >= 1) {
    stop("rank.tol must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.logical(center.tangent) || length(center.tangent) != 1L ||
      is.na(center.tangent)) {
    stop("center.tangent must be TRUE or FALSE.", call. = FALSE)
  }
  n <- length(riemobj$data)
  w <- riem_tangent_weights(weight, n)
  if (!is.null(weight) && any(weight > 0 & w == 0)) {
    stop("Positive weights underflowed during normalization; their relative scale is outside the numerical range.",
         call. = FALSE)
  }
  if (!center.tangent && sum(w > 0) < 2L) {
    stop("An uncentered fit requires at least two positive weights for its variance normalization.",
         call. = FALSE)
  }
  reference <- riem.mean(riemobj, weight = w, geometry = geom,
                         maxiter = maxiter, eps = eps)
  convergence <- reference$converged
  if (is.null(convergence) && !is.null(reference$diagnostics)) {
    convergence <- reference$diagnostics$converged
  }
  if (!isTRUE(convergence)) {
    stop("The reference mean did not converge; inspect riem.mean() before fitting tangent PCA.",
         call. = FALSE)
  }
  coordinates <- riem_tangent_spec(reference$mean, geom)
  z <- riem_tangent_encode(riemobj$data, coordinates)
  if (center.tangent) {
    # Center relative to an observed coordinate vector. This preserves exact
    # zero variation for repeated observations even when sum(w) rounds away
    # from one, without imposing a unit-dependent absolute rank threshold.
    coordinate.anchor <- z[which(w > 0)[1L], ]
    shifted.coordinates <- sweep(z, 2L, coordinate.anchor, "-")
    shifted.offset <- as.numeric(crossprod(w, shifted.coordinates))
    offset <- coordinate.anchor + shifted.offset
    centered <- sweep(shifted.coordinates, 2L, shifted.offset, "-")
  } else {
    coordinate.anchor <- rep(0, ncol(z))
    shifted.offset <- rep(0, ncol(z))
    offset <- rep(0, ncol(z))
    centered <- z
  }
  decomposition <- svd(sweep(centered, 1L, sqrt(w), "*"), nu = 0L)
  numerical.rank <- if (!length(decomposition$d) || decomposition$d[1L] == 0) {
    0L
  } else {
    sum(decomposition$d > rank.tol * decomposition$d[1L])
  }
  maximum.rank <- if (center.tangent) max(0L, sum(w > 0) - 1L) else sum(w > 0)
  numerical.rank <- min(numerical.rank, maximum.rank, geom$intrinsic_dimension)
  effective <- min(ndim, numerical.rank)
  if (effective < ndim) {
    warning(sprintf("Requested %s components, but numerical rank is %s; returning %s components.",
                    ndim, numerical.rank, effective), call. = FALSE)
  }
  directions <- decomposition$v[, seq_len(effective), drop = FALSE]
  # Fix signs for repeatable printing; repeated-eigenvalue subspaces remain nonunique.
  if (effective) {
    for (j in seq_len(effective)) {
      pivot <- which.max(abs(directions[, j]))
      if (directions[pivot, j] < 0) directions[, j] <- -directions[, j]
    }
  }
  scores <- centered %*% directions
  component.names <- if (effective > 0L) paste0("PC", seq_len(effective)) else character(0)
  colnames(scores) <- colnames(directions) <- component.names
  # 1 - sum(w^2) = 2 * sum_{i<j} w_i w_j for normalized weights.
  # Positive prefix products avoid subtracting two numbers close to one.
  normalization <- if (n > 1L) {
    2 * sum(w[-1L] * utils::head(cumsum(w), -1L))
  } else 0
  if (sum(w > 0) > 1L && (!is.finite(normalization) || normalization <= 0)) {
    stop("The weighted variance normalization is outside the positive numerical range.",
         call. = FALSE)
  }
  # Divide before squaring, so a small positive singular value is not lost to
  # underflow even when its variance relative to the denominator is ordinary.
  all.variance <- if (normalization > 0) {
    (decomposition$d / sqrt(normalization))^2
  } else rep(0, length(decomposition$d))
  if (any(!is.finite(all.variance))) {
    stop("Tangent variances are outside the finite numerical range.", call. = FALSE)
  }
  total.variance <- sum(all.variance)
  template <- riemobj
  template$data <- NULL
  output <- list(call = match.call(), center = reference$mean, embed = scores,
                 geometry = geom, coordinates = coordinates, offset = offset,
                 coordinate.anchor = coordinate.anchor, relative.offset = shifted.offset,
                 loadings = directions, singular.values = decomposition$d,
                 variance = all.variance[seq_len(effective)],
                 total.variance = total.variance,
                 explained.variance = if (total.variance > 0) {
                   all.variance[seq_len(effective)] / total.variance
                 } else rep(0, effective),
                 rank = numerical.rank, requested.dimension = ndim,
                 effective.dimension = effective, weight = w,
                 center.tangent = center.tangent, rank.tol = rank.tol,
                 input_template = template, schema_version = 1L,
                 package_version = as.character(utils::packageVersion("Riemann")),
                 diagnostics = list(converged = TRUE, mean = reference,
                                    approximation = "tangent PCA"))
  structure(output, class = c("riem_pga", "riemfit"))
}

#' @rdname riem.pga
#' @param object a fitted tangent PCA object.
#' @param newdata a compatible \code{riemdata} object.
#' @param ... additional arguments; currently unused.
#' @method predict riem_pga
#' @export
predict.riem_pga <- function(object, newdata, ...) {
  riem_tangent_check_fit(object)
  riem_check_newdata(object$input_template, newdata)
  z <- riem_tangent_encode(newdata$data, object$coordinates)
  shifted <- sweep(z, 2L, object$coordinate.anchor, "-")
  scores <- sweep(shifted, 2L, object$relative.offset, "-") %*% object$loadings
  dimnames(scores) <- list(NULL, colnames(object$embed))
  scores
}

#' Reconstruct Manifold Observations from Fitted Coordinates
#'
#' Reconstruct observations using the trained coordinate system of a model.
#' @param object a fitted model with a reconstruction method.
#' @param ... arguments passed to a reconstruction method.
#' @return a \code{riemdata} object in the fitted geometry and representation.
#' @export
riem.reconstruct <- function(object, ...) UseMethod("riem.reconstruct")

#' @rdname riem.pga
#' @param scores numeric matrix with one column per retained component; by
#'   default, the fitted training scores are reconstructed.
#' @method riem.reconstruct riem_pga
#' @export
riem.reconstruct.riem_pga <- function(object, scores = object$embed, ...) {
  riem_tangent_check_fit(object)
  if (is.numeric(scores) && is.null(dim(scores))) {
    if (object$effective.dimension == 1L) {
      scores <- matrix(scores, ncol = 1L)
    } else scores <- matrix(scores, nrow = 1L)
  }
  if (!is.matrix(scores) || !is.numeric(scores) || is.complex(scores) ||
      any(!is.finite(scores)) || ncol(scores) != object$effective.dimension ||
      nrow(scores) < 1L) {
    stop("scores must be a finite numeric matrix with one column per retained component.",
         call. = FALSE)
  }
  z <- sweep(scores %*% t(object$loadings), 2L, object$relative.offset, "+")
  z <- sweep(z, 2L, object$coordinate.anchor, "+")
  result <- object$input_template
  result$data <- riem_tangent_decode(z, object$coordinates)
  result$geometry <- object$geometry
  riem_validate_data(result)
  result
}

#' @rdname riem.pga
#' @method print riem_pga
#' @export
print.riem_pga <- function(x, ...) {
  cat("Tangent PCA on", x$geometry$manifold_id, "(", x$geometry$geometry_id, ")\n")
  cat(nrow(x$embed), "observations;", x$effective.dimension,
      "components retained; numerical rank", x$rank, "\n")
  invisible(x)
}

#' @rdname riem.pga
#' @param x a fitted tangent PCA object.
#' @method summary riem_pga
#' @export
summary.riem_pga <- function(object, ...) {
  riem_tangent_check_fit(object)
  importance <- rbind(variance = object$variance,
                      proportion = object$explained.variance,
                      cumulative = cumsum(object$explained.variance))
  colnames(importance) <- colnames(object$embed)
  structure(list(geometry = object$geometry, importance = importance,
                 rank = object$rank, requested.dimension = object$requested.dimension,
                 effective.dimension = object$effective.dimension,
                 approximation = "Variance explained in the fitted tangent representation"),
            class = "summary.riem_pga")
}

#' @rdname riem.pga
#' @method print summary.riem_pga
#' @export
print.summary.riem_pga <- function(x, ...) {
  cat(x$approximation, "\n")
  print(x$importance)
  invisible(x)
}

#' @rdname riem.pga
#' @param components two retained component indices to display.
#' @method plot riem_pga
#' @export
plot.riem_pga <- function(x, components = c(1L, 2L), ...) {
  riem_tangent_check_fit(x)
  if (x$effective.dimension == 0L) stop("No nonzero components to plot.", call. = FALSE)
  if (x$effective.dimension == 1L) {
    graphics::plot(seq_len(nrow(x$embed)), x$embed[, 1L],
                   xlab = "Observation", ylab = "PC1", ...)
  } else {
    if (!is.numeric(components) || length(components) != 2L ||
        any(!is.finite(components)) || any(components != floor(components)) ||
        any(components < 1L | components > x$effective.dimension)) {
      stop("components must identify two retained components.", call. = FALSE)
    }
    graphics::plot(x$embed[, components[1L]], x$embed[, components[2L]],
                   xlab = paste0("PC", components[1L]),
                   ylab = paste0("PC", components[2L]), ...)
  }
  invisible(x)
}
