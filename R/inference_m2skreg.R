#' Manifold-to-Scalar Kernel Regression
#'
#' Fits the Nadaraya--Watson smoother for manifold-valued predictors and finite
#' real scalar responses. For bandwidth \eqn{h>0}, weights at \eqn{x} are
#' proportional to \eqn{\exp\{-d(x,X_i)^2/(2h^2)\}}. The selected distance and
#' training observations are retained for prediction. Training fitted values
#' include the observation's own weight; they are not cross-validated predictions.
#'
#' @param riemobj A \code{riemdata} object containing the training predictors.
#' @param y A finite numeric response vector with one entry per predictor.
#' @param bandwidth One finite, strictly positive bandwidth.
#' @param geometry Geometry name or specification. \code{NULL} resolves the
#'   input's geometry, defaulting to its intrinsic geometry when unspecified.
#'   Legacy \code{"intrinsic"} and \code{"extrinsic"} aliases are accepted.
#'
#' @details Weights are normalized after subtracting the smallest squared
#'   distance in the exponent. The difference is evaluated without first
#'   squaring the distances, avoiding all-zero weights for small bandwidths or
#'   distant predictions. Numerically negligible relative weights can still
#'   underflow to zero. This is floating-point evaluation of the Gaussian
#'   smoother, not a change to a nearest-neighbor model.
#'
#'   The effective weight count and nearest training distance divided by the
#'   bandwidth describe the weights and distance scale. They are not confidence
#'   intervals or assurances of adequate statistical support.
#'
#' @return An object of class \code{m2skreg}. Legacy fields \code{ypred},
#'   \code{bandwidth}, and \code{inputs} are retained. Additional fields include
#'   resolved \code{geometry}, \code{call}, input dimensions, schema and package
#'   versions, and \code{training_diagnostics}. No distance matrix is retained.
#'
#' @examples
#' theta <- seq(0, pi / 2, length.out = 8)
#' X <- wrap.sphere(cbind(cos(theta), sin(theta)))
#' fit <- riem.m2skreg(X, sin(2 * theta), bandwidth = 0.3)
#' predict(fit, X)
#' summary(fit)
#'
#' @concept inference
#' @export
riem.m2skreg <- function(riemobj, y, bandwidth = 0.5, geometry = NULL) {
  riem_validate_data(riemobj)
  y <- riem_regression_response(y, length(riemobj$data))
  bandwidth <- riem_regression_bandwidth(bandwidth)
  geometry <- riem_resolve_geometry(riemobj, geometry, capability = "distance")
  distances <- basic_pdist(riemobj$name, riemobj$data, geometry$backend)
  result <- riem_kernel_predict_distances(distances, y, bandwidth)
  riem_regression_object(riemobj, y, bandwidth, geometry, result, match.call())
}

#' Prediction for Manifold-to-Scalar Kernel Regression
#'
#' Predicts under the fitted geometry and bandwidth without refitting or using
#' prediction-batch statistics. Training predictors and responses are required.
#'
#' @param object A fitted \code{m2skreg} object.
#' @param newdata A compatible \code{riemdata} object containing new predictors.
#' @param geometry Usually \code{NULL}, which uses the fitted geometry. An
#'   explicit value must resolve to the same geometry. An old serialized object
#'   without geometry metadata requires an explicit value and produces a
#'   warning, because its original geometry cannot be inferred reliably.
#' @param diagnostics Whether to return support diagnostics with predictions.
#' @param block_size Positive integer limiting the number of new observations
#'   in each cross-distance calculation. It does not change the fitted model.
#' @param ... Reserved for future arguments; unknown arguments are rejected.
#'
#' @return By default, a numeric vector with one prediction per observation.
#'   With \code{diagnostics=TRUE}, a list containing \code{prediction} and a data
#'   frame \code{diagnostics} with \code{effective_n}, \code{nearest_distance},
#'   and \code{nearest_over_bandwidth}. The last quantity may be infinite when
#'   the ratio is not representable even though the prediction remains finite.
#' @seealso \code{\link{riem.m2skreg}}
#' @concept inference
#' @method predict m2skreg
#' @export
predict.m2skreg <- function(object, newdata, geometry = NULL,
                           diagnostics = FALSE, block_size = 256L, ...) {
  if (length(list(...))) stop("Unknown prediction arguments.", call. = FALSE)
  riem_regression_validate_object(object)
  training <- object$inputs[[1L]]
  riem_check_newdata(training, newdata)
  y <- riem_regression_response(object$inputs[[2L]], length(training$data))
  bandwidth <- riem_regression_bandwidth(object$bandwidth)
  resolved <- riem_regression_prediction_geometry(object, geometry)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("'diagnostics' must be TRUE or FALSE.", call. = FALSE)
  }
  block_size <- riem_regression_integer(block_size, "block_size", 1L)
  nnew <- length(newdata$data)
  prediction <- numeric(nnew)
  support <- matrix(NA_real_, nnew, 3L,
                    dimnames = list(NULL, c("effective_n", "nearest_distance",
                                            "nearest_over_bandwidth")))
  for (first in seq.int(1L, nnew, by = block_size)) {
    ids <- seq.int(first, min(nnew, as.double(first) + block_size - 1))
    distances <- basic_pdist2(training$name, training$data, newdata$data[ids],
                             resolved$backend)
    result <- riem_kernel_predict_distances(distances, y, bandwidth)
    prediction[ids] <- result$prediction
    support[ids, ] <- as.matrix(result$diagnostics)
  }
  if (diagnostics) {
    return(list(prediction = prediction, diagnostics = as.data.frame(support)))
  }
  prediction
}

#' Methods for Scalar-Response Kernel Regression Fits
#'
#' @param object,x A fitted \code{m2skreg} object, or its summary for the summary
#'   printing method.
#' @param ... Additional arguments; currently unused.
#' @return \code{fitted} and \code{residuals} return numeric vectors; residuals
#'   are observed responses minus training fitted values. \code{summary}
#'   returns a \code{summary.m2skreg} list. Printing returns its input invisibly.
#' @name m2skreg-methods
#' @method fitted m2skreg
#' @export
fitted.m2skreg <- function(object, ...) {
  riem_regression_validate_object(object)
  object$ypred
}

#' @rdname m2skreg-methods
#' @method residuals m2skreg
#' @export
residuals.m2skreg <- function(object, ...) {
  riem_regression_validate_object(object)
  object$inputs[[2L]] - object$ypred
}

#' @rdname m2skreg-methods
#' @method summary m2skreg
#' @export
summary.m2skreg <- function(object, ...) {
  riem_regression_validate_object(object)
  residual <- stats::residuals(object)
  scale <- max(abs(residual))
  rmse <- if (scale == 0) 0 else if (!is.finite(scale)) Inf else
    sqrt(mean((residual / scale)^2)) * scale
  output <- list(call = object$call, nobs = length(object$inputs[[2L]]),
                 bandwidth = object$bandwidth, geometry = object$geometry,
                 training_rmse = rmse,
                 training_diagnostics = object$training_diagnostics,
                 errors = object$errors, fold_errors = object$fold_errors,
                 candidate_status = object$candidate_status,
                 cv_sse = if (!is.null(object$cv)) object$cv$selected_sse else NULL)
  structure(output, class = "summary.m2skreg")
}

#' @rdname m2skreg-methods
#' @method print m2skreg
#' @export
print.m2skreg <- function(x, ...) {
  riem_regression_validate_object(x)
  cat("Manifold-to-scalar Gaussian kernel regression\n")
  cat("Observations:", length(x$inputs[[2L]]), "  Bandwidth:", x$bandwidth, "\n")
  cat("Geometry:", if (is.null(x$geometry)) "unknown (legacy object)" else
    x$geometry$geometry_id, "\n")
  if (!is.null(x$cv)) cat("Cross-validation:", length(unique(x$foldid)),
                         "folds; selected SSE:", x$cv$selected_sse, "\n")
  invisible(x)
}

#' @rdname m2skreg-methods
#' @method print summary.m2skreg
#' @export
print.summary.m2skreg <- function(x, ...) {
  cat("Manifold-to-scalar Gaussian kernel regression\n")
  cat("Observations:", x$nobs, "  Bandwidth:", x$bandwidth, "\n")
  cat("Geometry:", if (is.null(x$geometry)) "unknown (legacy object)" else
    x$geometry$geometry_id, "\n")
  cat("Training RMSE (includes self-weights):", x$training_rmse, "\n")
  if (!is.null(x$cv_sse)) cat("Selected cross-validation SSE:", x$cv_sse, "\n")
  invisible(x)
}
