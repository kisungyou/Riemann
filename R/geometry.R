#' Inspect and Resolve a Statistical Geometry
#'
#' Geometry is resolved once when a model is fitted and retained for prediction.
#' For SPD matrices, \code{"intrinsic"} means \code{"affine_invariant"} and
#' \code{"extrinsic"} means \code{"log_euclidean"}. On the sphere the corresponding
#' names are \code{"round"} and \code{"chordal"}. Euclidean aliases coincide.
#' Landmark geometry currently identifies configurations under the full orthogonal
#' group, including reflections. A capability records an available operation;
#' its validation status must also be considered before interpreting results.
#'
#' @param riemobj A wrapped \code{riemdata} object.
#' @param geometry A geometry name, a saved \code{riem_geometry} specification,
#'   or \code{NULL} for the default intrinsic geometry.
#' @return \code{riem.geometry} returns a serializable geometry specification.
#'   \code{riem.capabilities} returns the operation registry as a data frame.
#' @examples
#' x <- wrap.spd(list(diag(2), 2 * diag(2)))
#' riem.geometry(x, "log_euclidean")
#' riem.capabilities()
#' @export
riem.geometry <- function(riemobj, geometry = NULL) {
  riem_resolve_geometry(riemobj, geometry)
}

#' @rdname riem.geometry
#' @export
riem.capabilities <- function() {
  manifolds <- c("spd", "sphere", "euclidean", "landmark", "grassmann",
                 "stiefel", "rotation", "multinomial", "spdk", "correlation")
  registry <- expand.grid(manifold_id = manifolds,
                          backend = c("intrinsic", "extrinsic"),
                          stringsAsFactors = FALSE)
  registry$geometry_id <- mapply(riem_geometry_id, registry$manifold_id,
                                 registry$backend, USE.NAMES = FALSE)
  registry$distance <- TRUE
  unsupported <- registry$manifold_id == "correlation" |
    (registry$manifold_id == "stiefel" & registry$backend == "intrinsic") |
    (registry$manifold_id == "spdk" & registry$backend == "extrinsic")
  registry$distance[unsupported] <- FALSE
  registry$mean <- !(registry$manifold_id == "landmark" &
                       registry$backend == "extrinsic")
  registry$median <- registry$mean
  registry$mean[unsupported] <- FALSE
  registry$median[unsupported] <- FALSE
  registry$tangent <- registry$manifold_id %in% c("spd", "euclidean") |
    (registry$manifold_id %in% c("sphere", "landmark") &
       registry$backend == "intrinsic")
  registry$validation_status <- ifelse(
    registry$manifold_id %in% c("spd", "sphere", "euclidean"),
    "core_contract", "audit_pending")
  registry$validation_status[unsupported] <- "restricted_inconsistent_geometry"
  registry$median_estimand <- ifelse(
    registry$backend == "extrinsic" &
      !registry$manifold_id %in% c("spd", "euclidean"),
    "projected_ambient_median", "frechet_median")
  registry
}

riem_geometry_id <- function(manifold, backend) {
  if (manifold == "spd") {
    return(if (backend == "intrinsic") "affine_invariant" else "log_euclidean")
  }
  if (manifold == "sphere") {
    return(if (backend == "intrinsic") "round" else "chordal")
  }
  if (manifold == "euclidean") return("euclidean")
  if (manifold == "landmark") {
    return(if (backend == "intrinsic") "shape_orthogonal" else "procrustes_chordal")
  }
  definitions <- list(
    grassmann = c(intrinsic = "principal_angles", extrinsic = "projector_chordal"),
    stiefel = c(intrinsic = "stiefel_intrinsic_unavailable", extrinsic = "frame_chordal"),
    rotation = c(intrinsic = "rotation_frobenius", extrinsic = "rotation_chordal"),
    multinomial = c(intrinsic = "fisher_rao", extrinsic = "sqrt_chordal"),
    spdk = c(intrinsic = "factor_procrustes", extrinsic = "factor_extrinsic_unavailable"),
    correlation = c(intrinsic = "correlation_quotient_unavailable",
                    extrinsic = "correlation_extrinsic_unavailable"))
  unname(definitions[[manifold]][backend])
}

riem_resolve_geometry <- function(riemobj, geometry = NULL, capability = NULL) {
  riem_validate_data(riemobj)
  if (is.null(geometry) && inherits(riemobj$geometry, "riem_geometry")) geometry <- riemobj$geometry
  saved <- inherits(geometry, "riem_geometry")
  if (saved) {
    if (!identical(geometry$schema_version, 1L) ||
        !identical(geometry$manifold_id, riemobj$name) ||
        !identical(as.integer(geometry$representation_dim), as.integer(riemobj$size))) {
      stop("Saved geometry has an incompatible schema, manifold, or dimensions.", call. = FALSE)
    }
    selector <- geometry$geometry_id
  } else {
    selector <- if (is.null(geometry)) "intrinsic" else geometry
  }
  if (!is.character(selector) || length(selector) != 1L || is.na(selector)) {
    stop("geometry must be one name or a saved riem_geometry specification.", call. = FALSE)
  }
  selector <- tolower(selector)
  # Accept names used by early development snapshots as compatibility aliases.
  if (selector %in% c("native_intrinsic", "native_extrinsic") && !saved) {
    selector <- sub("^native_", "", selector)
  }
  aliases <- c(airm = "affine_invariant", logeuclidean = "log_euclidean")
  if (selector %in% names(aliases)) selector <- unname(aliases[selector])
  registry <- riem.capabilities()
  rows <- registry[registry$manifold_id == riemobj$name, , drop = FALSE]
  if (selector %in% c("intrinsic", "extrinsic")) {
    rows <- rows[rows$backend == selector, , drop = FALSE]
  } else {
    rows <- rows[rows$geometry_id == selector, , drop = FALSE]
  }
  if (!nrow(rows)) stop("Unsupported geometry for ", riemobj$name, ": ", selector, call. = FALSE)
  row <- rows[1L, , drop = FALSE]
  if (riemobj$name == "euclidean") row$backend <- "intrinsic"
  if (!is.null(capability) &&
      (!capability %in% names(row) || !isTRUE(row[[capability]][1L]))) {
    stop("Geometry '", row$geometry_id, "' does not support ", capability,
         " with a compatible implementation.", call. = FALSE)
  }
  d <- as.integer(riemobj$size)
  q <- switch(riemobj$name, spd = d[1] * (d[1] + 1) / 2,
    sphere = prod(d) - 1, euclidean = prod(d),
    landmark = (d[1] - 1) * d[2] - 1 - d[2] * (d[2] - 1) / 2,
    grassmann = d[2] * (d[1] - d[2]),
    stiefel = d[1] * d[2] - d[2] * (d[2] + 1) / 2,
    rotation = d[1] * (d[1] - 1) / 2,
    multinomial = prod(d) - 1, correlation = d[1] * (d[1] - 1) / 2,
    spdk = d[1] * d[2] - d[2] * (d[2] - 1) / 2)
  parameters <- if (riemobj$name == "landmark") list(reflections = TRUE) else list()
  if (saved && !identical(parameters, geometry$parameters)) {
    stop("Unsupported or incompatible saved geometry parameters.", call. = FALSE)
  }
  result <- structure(list(manifold_id = riemobj$name, geometry_id = row$geometry_id,
    backend = row$backend, representation_dim = d, intrinsic_dimension = q,
    parameters = parameters, representation = "matrix",
    tangent_representation = if (row$geometry_id == "log_euclidean") "log_chart" else "tangent",
    distance_kind = if (row$validation_status == "core_contract") "metric" else "unclassified",
    distance_scale = "distance", capabilities = names(row)[vapply(row, isTRUE, logical(1))],
    numerical_policy = list(cut_locus = "error", data_repair = "none"),
    validation_status = row$validation_status, schema_version = 1L), class = "riem_geometry")
  result$tangent_representation <- switch(riemobj$name,
    rotation = "body_skew", grassmann = "horizontal_basis_lift",
    spdk = "horizontal_factor_lift", landmark = "horizontal_preshape_tangent",
    result$tangent_representation)
  if (riemobj$name == "spdk") result$representation <- "YYt_factor"
  if (saved) {
    fields <- c("manifold_id", "geometry_id", "backend", "representation_dim",
                "parameters", "representation", "tangent_representation",
                "distance_scale", "schema_version")
    if (!identical(geometry[fields], result[fields])) {
      stop("Saved geometry contains inconsistent or missing defining metadata; refit with a supported geometry.",
           call. = FALSE)
    }
  }
  result
}

#' @method print riem_geometry
#' @export
print.riem_geometry <- function(x, ...) {
  cat("Geometry:", x$geometry_id, "on", x$manifold_id, "\n")
  cat("Representation:", paste(x$representation_dim, collapse = " x "),
      "; intrinsic dimension:", x$intrinsic_dimension, "\n")
  cat("Validation:", x$validation_status, "\n")
  invisible(x)
}

riem_validate_data <- function(riemobj) {
  if (!inherits(riemobj, "riemdata") || !is.list(riemobj$data) ||
      !length(riemobj$data) || !is.character(riemobj$name) ||
      length(riemobj$name) != 1L || is.na(riemobj$name)) {
    stop("Expected a nonempty riemdata object created by a wrap function.", call. = FALSE)
  }
  dims <- riemobj$size
  if (!is.numeric(dims) || length(dims) != 2L || anyNA(dims) ||
      any(!is.finite(dims)) || any(dims < 1) || any(dims != trunc(dims))) {
    stop("Invalid riemdata representation dimensions.", call. = FALSE)
  }
  for (i in seq_along(riemobj$data)) {
    x <- riemobj$data[[i]]
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || any(!is.finite(x)) ||
        !identical(as.integer(dim(x)), as.integer(dims))) {
      stop("Observation ", i, " must be a finite real matrix of the recorded dimensions.", call. = FALSE)
    }
    if (!is.null(riemobj$dimnames) && !identical(dimnames(x), riemobj$dimnames)) {
      stop("Observation ", i, " has inconsistent feature or landmark ordering.", call. = FALSE)
    }
    if (!identical(dimnames(x), dimnames(riemobj$data[[1L]]))) {
      stop("Observation ", i, " has inconsistent feature or landmark ordering.", call. = FALSE)
    }
    if (riemobj$name == "spd") check_spd(x, i)
    if (riemobj$name == "correlation") check_corr(x, i)
    if (riemobj$name == "rotation") single_rotcheck(x, i)
    if (riemobj$name %in% c("grassmann", "stiefel") &&
        (nrow(x) < ncol(x) || norm(crossprod(x) - diag(ncol(x)), "F") > 1e-8)) {
      stop("Frame observations must have orthonormal columns; use the appropriate wrapper.", call. = FALSE)
    }
    if (riemobj$name == "spdk" && !riem_full_column_rank(x)) {
      stop("Fixed-rank factors must have full column rank.", call. = FALSE)
    }
    if (riemobj$name == "multinomial" &&
        (length(x) < 2L || any(x <= 0) || abs(sum(x) - 1) > 1e-10)) {
      stop("Simplex observations must be strictly positive and sum to one.", call. = FALSE)
    }
    if (riemobj$name == "landmark" &&
        (abs(sqrt(sum(x*x)) - 1) > 1e-8 || max(abs(colMeans(x))) > 1e-8 ||
         !riem_full_column_rank(x))) {
      stop("Landmark observations must be centered unit preshapes; use wrap.landmark().", call. = FALSE)
    }
    if (riemobj$name == "sphere" && abs(sqrt(sum(x * x)) - 1) > 1e-8) {
      stop("Sphere observations must have unit norm; use wrap.sphere().", call. = FALSE)
    }
  }
  invisible(TRUE)
}

riem_check_newdata <- function(training, newdata) {
  riem_validate_data(newdata)
  if (!identical(training$name, newdata$name) ||
      !identical(as.integer(training$size), as.integer(newdata$size))) {
    stop("newdata has an incompatible manifold or representation dimensions.", call. = FALSE)
  }
  for (field in c("dimnames", "feature_names", "landmark_names", "representation")) {
    if (!identical(training[[field]], newdata[[field]])) {
      stop("newdata has incompatible ", field, "; preserve training feature order.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

riem_input_contract <- function(x) {
  x$data <- NULL
  x
}

riem_vector_input <- function(input, manifold) {
  feature_names <- if (is.matrix(input)) colnames(input) else NULL
  if (is.matrix(input)) {
    if (!nrow(input) || !ncol(input)) stop("Input must be nonempty.", call. = FALSE)
    data <- lapply(seq_len(nrow(input)), function(i) input[i, ])
  } else if (is.list(input)) data <- input
  else stop("Input must be a matrix of row observations or a list of vectors.", call. = FALSE)
  if (!check_list_eqsize(data)) stop("Observations must be nonempty vectors of the same length.", call. = FALSE)
  result <- lapply(seq_along(data), function(i) {
    x <- data[[i]]
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || any(!is.finite(x))) {
      stop("Observation ", i, " must be a finite real vector.", call. = FALSE)
    }
    if (manifold == "sphere") {
      if (length(x) < 2L || max(abs(x)) == 0) {
        stop("Sphere vectors must have at least two entries and nonzero norm.", call. = FALSE)
      }
      x <- x / max(abs(x))
      x <- x / sqrt(sum(x * x))
    }
    nm <- if (is.null(feature_names)) names(x) else feature_names
    matrix(x, ncol = 1L, dimnames = if (is.null(nm)) NULL else list(nm, NULL))
  })
  if (!all(vapply(result, function(x) identical(dimnames(x), dimnames(result[[1L]])), logical(1)))) {
    stop("Observations must have consistent feature ordering, including whether names are supplied.", call. = FALSE)
  }
  result
}
