# Internal scalar-response regression utilities. Geometry calculations are
# delegated to the shared registry and native distance implementations.

riem_regression_response <- function(y, n) {
  if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) ||
      length(y) != n || any(!is.finite(y))) {
    stop("'y' must be a finite real numeric vector with one response per observation.",
         call. = FALSE)
  }
  as.double(y)
}

riem_regression_bandwidth <- function(bandwidth, multiple = FALSE) {
  if (!is.numeric(bandwidth) || is.complex(bandwidth) ||
      !is.null(dim(bandwidth)) || !length(bandwidth) ||
      (!multiple && length(bandwidth) != 1L) ||
      any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
    stop(if (multiple) "'bandwidths' must contain finite, strictly positive numbers." else
      "'bandwidth' must be one finite, strictly positive number.", call. = FALSE)
  }
  as.double(bandwidth)
}

riem_regression_integer <- function(value, name, lower, upper = .Machine$integer.max) {
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
      !is.finite(value) || value != floor(value) || value < lower || value > upper) {
    stop(sprintf("'%s' must be an integer between %s and %s.", name, lower, upper),
         call. = FALSE)
  }
  as.integer(value)
}

riem_kernel_weights <- function(distances, bandwidth) {
  if (!is.numeric(distances) || is.complex(distances) || !length(distances) ||
      any(!is.finite(distances)) || any(distances < 0)) {
    stop("Kernel regression requires finite, nonnegative distances.", call. = FALSE)
  }
  nearest <- min(distances)
  delta <- distances - nearest
  weights <- rep(1, length(distances))
  positive <- delta > 0
  # (d^2 - dmin^2)/(2h^2) = (d-dmin)(d+dmin)/(2h^2).
  # Compute its logarithm to avoid overflow in either squares or d + dmin.
  if (any(positive)) {
    d <- distances[positive]
    log_gap <- log(delta[positive]) + log(d) + log1p(nearest / d) -
      log(2) - 2 * log(bandwidth)
    weights[positive] <- exp(-exp(log_gap))
  }
  weights / sum(weights)
}

riem_kernel_predict_distances <- function(distances, y, bandwidth) {
  if (!is.matrix(distances) || nrow(distances) != length(y) ||
      !ncol(distances) || !nrow(distances)) {
    stop("Distance dimensions do not match the training responses and predictions.",
         call. = FALSE)
  }
  if (!is.numeric(distances) || is.complex(distances) ||
      any(!is.finite(distances)) || any(distances < 0)) {
    stop("Kernel regression requires finite, nonnegative distances.", call. = FALSE)
  }
  prediction <- numeric(ncol(distances))
  support <- matrix(0, ncol(distances), 3L,
                    dimnames = list(NULL, c("effective_n", "nearest_distance",
                                            "nearest_over_bandwidth")))
  response_scale <- max(abs(y))
  scaled_y <- if (response_scale == 0) y else y / response_scale
  for (j in seq_len(ncol(distances))) {
    weights <- riem_kernel_weights(distances[, j], bandwidth)
    # Bound roundoff by the convex hull; even maximal finite responses remain
    # finite when the exact weighted average equals an endpoint.
    average <- min(max(scaled_y), max(min(scaled_y), sum(weights * scaled_y)))
    prediction[j] <- average * response_scale
    nearest <- min(distances[, j])
    support[j, ] <- c(1 / sum(weights^2), nearest, nearest / bandwidth)
  }
  if (any(!is.finite(prediction))) {
    stop("The weighted response could not be represented as a finite number.", call. = FALSE)
  }
  list(prediction = prediction, diagnostics = as.data.frame(support))
}

riem_regression_object <- function(riemobj, y, bandwidth, geometry, result, call) {
  structure(list(ypred = result$prediction, bandwidth = bandwidth,
                 inputs = list(riemobj, y), call = call,
                 task = "scalar_response_kernel_regression", geometry = geometry,
                 input_dimensions = riemobj$size, nobs = length(y),
                 training_diagnostics = result$diagnostics,
                 controls = list(kernel = "gaussian", weight_normalization = "sum_to_one"),
                 converged = TRUE, termination = "closed_form", iterations = 0L,
                 package_version = as.character(utils::packageVersion("Riemann")),
                 schema_version = 1L), class = c("m2skreg", "riemfit"))
}

riem_regression_validate_object <- function(object) {
  if (!inherits(object, "m2skreg") || !is.list(object$inputs) ||
      length(object$inputs) != 2L || !is.numeric(object$ypred) || is.complex(object$ypred) ||
      length(object$ypred) != length(object$inputs[[2L]]) ||
      any(!is.finite(object$ypred))) {
    stop("Invalid 'm2skreg' object: finite fitted values and training inputs are required.",
         call. = FALSE)
  }
  riem_validate_data(object$inputs[[1L]])
  riem_regression_response(object$inputs[[2L]], length(object$inputs[[1L]]$data))
  riem_regression_bandwidth(object$bandwidth)
  if (!is.null(object$schema_version) &&
      !identical(object$schema_version, 1L)) {
    stop("Unsupported 'm2skreg' object schema; use a compatible Riemann version.", call. = FALSE)
  }
  invisible(TRUE)
}

riem_regression_prediction_geometry <- function(object, geometry) {
  training <- object$inputs[[1L]]
  if (is.null(object$geometry)) {
    if (is.null(geometry)) {
      stop("This legacy fit has no geometry metadata. Supply its original geometry explicitly or refit it.",
           call. = FALSE)
    }
    warning("Legacy fit: using the explicitly supplied geometry; the original geometry cannot be verified.",
            call. = FALSE)
    return(riem_resolve_geometry(training, geometry, capability = "distance"))
  }
  stored <- riem_resolve_geometry(training, object$geometry, capability = "distance")
  if (!is.null(geometry)) {
    proposed <- riem_resolve_geometry(training, geometry, capability = "distance")
    fields <- c("manifold_id", "geometry_id", "backend", "parameters")
    if (!identical(stored[fields], proposed[fields])) {
      stop("Prediction geometry must match the fitted geometry; refit to change it.", call. = FALSE)
    }
  }
  stored
}
