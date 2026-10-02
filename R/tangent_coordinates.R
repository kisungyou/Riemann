# Internal metric-preserving coordinates for tangent PCA.

riem_tangent_weights <- function(weight, n) {
  if (is.null(weight)) return(rep(1 / n, n))
  if (!is.numeric(weight) || is.complex(weight) || length(weight) != n ||
      any(!is.finite(weight)) || any(weight < 0) || !any(weight > 0)) {
    stop("weight must contain finite nonnegative values with positive total mass.", call. = FALSE)
  }
  weight <- weight / max(weight)
  weight / sum(weight)
}

riem_tangent_symmetric_function <- function(x, fun, positive = FALSE) {
  x <- x / 2 + t(x) / 2
  spectral <- eigen(x, symmetric = TRUE)
  if (positive && any(spectral$values <= 0)) {
    stop("A tangent-coordinate matrix is not positive definite.", call. = FALSE)
  }
  values <- fun(spectral$values)
  if (any(!is.finite(values))) {
    stop("The requested matrix function is outside the finite numerical range.", call. = FALSE)
  }
  output <- tcrossprod(sweep(spectral$vectors, 2L, values, "*"), spectral$vectors)
  output / 2 + t(output) / 2
}

riem_tangent_svec <- function(x) {
  keep <- upper.tri(x, diag = TRUE)
  x[keep] * ifelse(row(x)[keep] == col(x)[keep], 1, sqrt(2))
}

riem_tangent_unsvec <- function(x, p) {
  output <- matrix(0, p, p)
  keep <- upper.tri(output, diag = TRUE)
  output[keep] <- x / ifelse(row(output)[keep] == col(output)[keep], 1, sqrt(2))
  output + t(output) - diag(diag(output), nrow = p, ncol = p)
}

riem_tangent_sphere_log <- function(center, x) {
  center <- as.numeric(center)
  x <- as.numeric(x)
  # The half-angle formula retains tiny angles when dot products round to one.
  angle <- 2 * atan2(sqrt(sum((center - x)^2)), sqrt(sum((center + x)^2)))
  if (sqrt(sum((center + x)^2)) <= 64 * .Machine$double.eps) {
    stop("The tangent logarithm is nonunique at an antipode.", call. = FALSE)
  }
  if (angle == 0) return(rep(0, length(center)))
  direction <- x - sum(center * x) * center
  direction <- direction - sum(direction * center) * center
  length.direction <- sqrt(sum(direction^2))
  if (length.direction == 0) {
    stop("The tangent logarithm cannot resolve this direction numerically.", call. = FALSE)
  }
  direction * (angle / length.direction)
}

riem_tangent_sphere_exp <- function(center, tangent, limit = pi) {
  center <- as.numeric(center)
  tangent <- as.numeric(tangent)
  theta <- sqrt(sum(tangent^2))
  if (!is.finite(theta) || theta >= limit) {
    stop("Reconstruction exceeds the supported injectivity neighborhood.", call. = FALSE)
  }
  if (theta == 0) return(center)
  result <- cos(theta) * center + (sin(theta) / theta) * tangent
  result / sqrt(sum(result^2))
}

riem_tangent_landmark_align <- function(center, x, tol = sqrt(.Machine$double.eps)) {
  decomposition <- svd(crossprod(center, x))
  # At orthogonal spans the whole cross-Gram matrix is roundoff-sized. Its
  # condition number alone cannot distinguish this cut locus from a valid lift.
  scale <- norm(center, "2") * norm(x, "2")
  if (min(decomposition$d) <= tol * scale) {
    stop("Landmark alignment is singular or numerically nonunique at this reference.",
         call. = FALSE)
  }
  x %*% decomposition$v %*% t(decomposition$u)
}

riem_tangent_landmark_regular <- function(x) {
  values <- svd(x, nu = 0L, nv = 0L)$d
  length(values) == ncol(x) && min(values) > sqrt(.Machine$double.eps) * max(values)
}

riem_tangent_spec <- function(center, geometry) {
  spec <- list(center = center, geometry = geometry, schema_version = 1L,
               dimensions = dim(center))
  id <- geometry$geometry_id
  manifold <- geometry$manifold_id
  if (manifold == "spd" && id == "affine_invariant") {
    spec$type <- "spd_affine_svec"
    spec$root <- riem_tangent_symmetric_function(center, sqrt, positive = TRUE)
    spec$inverse.root <- riem_tangent_symmetric_function(center, function(x) 1 / sqrt(x),
                                                        positive = TRUE)
    spec$coordinate.dimension <- nrow(center) * (nrow(center) + 1L) / 2L
  } else if (manifold == "spd" && id == "log_euclidean") {
    spec$type <- "spd_log_svec"
    spec$log.center <- riem_tangent_symmetric_function(center, log, positive = TRUE)
    spec$coordinate.dimension <- nrow(center) * (nrow(center) + 1L) / 2L
  } else if (manifold == "euclidean" && id == "euclidean") {
    spec$type <- "euclidean"
    spec$coordinate.dimension <- length(center)
  } else if (manifold == "sphere" && id == "round") {
    spec$type <- "sphere_orthonormal"
    # A Householder reflection maps the first coordinate vector to center.
    unit <- as.numeric(center)
    v <- unit
    sign.first <- if (v[1L] >= 0) 1 else -1
    v[1L] <- v[1L] + sign.first
    h <- diag(length(unit)) - 2 * tcrossprod(v) / sum(v^2)
    spec$basis <- h[, -1L, drop = FALSE]
    spec$coordinate.dimension <- ncol(spec$basis)
  } else if (manifold == "landmark" && id == "shape_orthogonal") {
    if (!riem_tangent_landmark_regular(center)) {
      stop("Landmark tangent PCA requires a regular full-column-rank reference shape.",
           call. = FALSE)
    }
    spec$type <- "landmark_horizontal"
    spec$alignment.group <- "O(p)"
    spec$coordinate.dimension <- length(center)
    spec$intrinsic.dimension <- (nrow(center) - 1L) * ncol(center) - 1L -
      ncol(center) * (ncol(center) - 1L) / 2L
  } else {
    stop("Tangent PCA has no validated coordinates for this manifold and geometry.",
         call. = FALSE)
  }
  spec
}

riem_tangent_encode <- function(data, spec) {
  output <- matrix(0, length(data), spec$coordinate.dimension)
  for (i in seq_along(data)) {
    x <- data[[i]]
    value <- switch(spec$type,
      spd_affine_svec = {
        relative <- spec$inverse.root %*% x %*% spec$inverse.root
        riem_tangent_svec(riem_tangent_symmetric_function(relative, log, positive = TRUE))
      },
      spd_log_svec = {
        riem_tangent_svec(riem_tangent_symmetric_function(x, log, positive = TRUE) -
                           spec$log.center)
      },
      euclidean = as.numeric(x - spec$center),
      sphere_orthonormal = as.numeric(crossprod(spec$basis,
                                               riem_tangent_sphere_log(spec$center, x))),
      landmark_horizontal = {
        aligned <- riem_tangent_landmark_align(spec$center, x)
        riem_tangent_sphere_log(spec$center, aligned)
      },
      stop("Unsupported stored coordinate convention.", call. = FALSE))
    output[i, ] <- value
  }
  if (any(!is.finite(output))) stop("Nonfinite tangent coordinates.", call. = FALSE)
  output
}

riem_tangent_decode <- function(coordinates, spec) {
  lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- coordinates[i, ]
    result <- switch(spec$type,
      spd_affine_svec = {
        tangent <- riem_tangent_unsvec(coordinate, nrow(spec$center))
        spec$root %*% riem_tangent_symmetric_function(tangent, exp) %*% spec$root
      },
      spd_log_svec = {
        tangent <- riem_tangent_unsvec(coordinate, nrow(spec$center))
        riem_tangent_symmetric_function(spec$log.center + tangent, exp)
      },
      euclidean = spec$center + matrix(coordinate, nrow = spec$dimensions[1L],
                                       ncol = spec$dimensions[2L]),
      sphere_orthonormal = matrix(riem_tangent_sphere_exp(spec$center,
                                                        spec$basis %*% coordinate),
                                  nrow = spec$dimensions[1L], ncol = spec$dimensions[2L]),
      landmark_horizontal = {
        shape <- matrix(riem_tangent_sphere_exp(spec$center, coordinate, limit = pi / 2),
                        nrow = spec$dimensions[1L], ncol = spec$dimensions[2L])
        shape <- sweep(shape, 2L, colMeans(shape), "-")
        shape <- shape / sqrt(sum(shape^2))
        if (!riem_tangent_landmark_regular(shape)) {
          stop("Reconstruction reaches an unsupported singular landmark shape.", call. = FALSE)
        }
        # Check that this is still in the trained local quotient chart.
        aligned <- riem_tangent_landmark_align(spec$center, shape)
        if (sqrt(sum((aligned - shape)^2)) > 1e-7) {
          stop("Reconstruction leaves the trained landmark alignment neighborhood.", call. = FALSE)
        }
        shape
      },
      stop("Unsupported stored coordinate convention.", call. = FALSE))
    dimnames(result) <- dimnames(spec$center)
    if (any(!is.finite(result))) stop("Nonfinite reconstruction.", call. = FALSE)
    result
  })
}

riem_tangent_check_fit <- function(object) {
  if (!inherits(object, "riem_pga") || !identical(object$schema_version, 1L) ||
      !is.list(object$coordinates) || !identical(object$coordinates$schema_version, 1L) ||
      is.null(object$geometry) || is.null(object$loadings) || is.null(object$offset) ||
      is.null(object$coordinate.anchor) || is.null(object$relative.offset) ||
      is.null(object$input_template)) {
    stop("This is not a supported fitted tangent PCA object; refit with riem.pga().",
         call. = FALSE)
  }
  template <- object$input_template
  template$data <- list(object$coordinates$center)
  geom <- riem_resolve_geometry(template, object$geometry, capability = "tangent")
  stored <- riem_resolve_geometry(template, object$coordinates$geometry, capability = "tangent")
  if (!identical(geom, stored) ||
      !is.matrix(object$loadings) || any(!is.finite(object$loadings)) ||
      nrow(object$loadings) != object$coordinates$coordinate.dimension ||
      length(object$offset) != nrow(object$loadings) || any(!is.finite(object$offset)) ||
      length(object$coordinate.anchor) != nrow(object$loadings) ||
      any(!is.finite(object$coordinate.anchor)) ||
      length(object$relative.offset) != nrow(object$loadings) ||
      any(!is.finite(object$relative.offset)) ||
      !identical(object$offset, object$coordinate.anchor + object$relative.offset)) {
    stop("Stored tangent PCA coordinates or geometry are incompatible; refit the model.", call. = FALSE)
  }
  invisible(TRUE)
}
