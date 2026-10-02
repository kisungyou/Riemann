riem_kmeans_integer <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > maximum) {
    stop(name, " must be an integer in [", minimum, ", ", maximum, "].", call. = FALSE)
  }
  as.integer(x)
}

riem_kmeans_distances <- function(data, centers, geometry) {
  distances <- basic_pdist2(geometry$manifold_id, data, centers, geometry$backend)
  distances <- matrix(distances, nrow = length(data), ncol = length(centers))
  if (any(!is.finite(distances)) || any(distances < 0)) {
    stop("Clustering encountered invalid or nonfinite distances.", call. = FALSE)
  }
  distances
}

riem_kmeans_assign <- function(distances) {
  max.col(-distances, ties.method = "first")
}

riem_kmeans_objective <- function(distances, labels) {
  objective <- sum(distances[cbind(seq_len(nrow(distances)), labels)]^2)
  if (!is.finite(objective)) stop("Clustering objective overflowed.", call. = FALSE)
  objective
}

riem_kmeans_initialize <- function(data, k, init, geometry) {
  n <- length(data)
  if (is.numeric(init)) return(as.integer(init))
  if (init == "random") return(sample.int(n, k))
  selected <- integer(k)
  selected[1L] <- sample.int(n, 1L)
  nearest <- riem_kmeans_distances(data, data[selected[1L]], geometry)[, 1L]
  if (k > 1L) {
    for (j in 2:k) {
      nearest[selected[seq_len(j - 1L)]] <- 0
      scale <- max(nearest)
      if (scale == 0) stop("k exceeds the number of distinct observations in this geometry.",
                            call. = FALSE)
      # Squared distances, scaled before squaring to avoid overflow.
      probability <- (nearest / scale)^2
      selected[j] <- sample.int(n, 1L, prob = probability)
      nearest <- pmin(nearest, riem_kmeans_distances(data, data[selected[j]], geometry)[, 1L])
    }
  }
  selected
}

riem_kmeans_repair_empty <- function(labels, distances, k) {
  counts <- tabulate(labels, nbins = k)
  empty <- which(counts == 0L)
  if (!length(empty)) return(labels)
  residual <- distances[cbind(seq_along(labels), labels)]
  for (cluster in empty) {
    eligible <- which(counts[labels] > 1L & residual > 0)
    if (!length(eligible)) {
      stop("Empty clusters cannot be separated: fewer distinct locations than requested clusters.",
           call. = FALSE)
    }
    # Stable tie rule: the lowest observation index among maximal residuals.
    chosen <- eligible[which.max(residual[eligible])]
    counts[labels[chosen]] <- counts[labels[chosen]] - 1L
    labels[chosen] <- cluster
    counts[cluster] <- 1L
    residual[chosen] <- 0
  }
  labels
}

riem_kmeans_mean <- function(riemobj, indices, geometry, mean.maxiter, mean.eps) {
  subset <- riemobj
  subset$data <- riemobj$data[indices]
  fit <- riem.mean(subset, geometry = geometry, maxiter = mean.maxiter, eps = mean.eps)
  if (!isTRUE(fit$converged)) {
    stop("A cluster mean failed to converge (", fit$termination, ").", call. = FALSE)
  }
  fit$mean
}

riem_kmeans_start <- function(riemobj, k, geometry, init, algorithm,
                              maxiter, eps, mean.maxiter, mean.eps) {
  seeds <- riem_kmeans_initialize(riemobj$data, k, init, geometry)
  centers <- riemobj$data[seeds]
  distances <- riem_kmeans_distances(riemobj$data, centers, geometry)
  labels <- riem_kmeans_assign(distances)
  objective <- riem_kmeans_objective(distances, labels)
  history <- objective
  converged <- FALSE
  termination <- "iteration_limit"
  empty.repairs <- 0L
  for (iteration in seq_len(maxiter)) {
    old.centers <- centers
    old.distances <- distances
    old.labels <- labels
    old.objective <- objective
    empty.repairs <- empty.repairs + sum(tabulate(labels, k) == 0L)
    labels <- riem_kmeans_repair_empty(labels, distances, k)
    centers <- lapply(seq_len(k), function(cluster) {
      riem_kmeans_mean(riemobj, which(labels == cluster), geometry, mean.maxiter, mean.eps)
    })
    if (algorithm == "macqueen") {
      # Sequential reassignment with exact/iterative Frechet updates of affected clusters.
      counts <- tabulate(labels, k)
      for (i in seq_along(labels)) {
        point.distances <- riem_kmeans_distances(riemobj$data[i], centers, geometry)
        destination <- riem_kmeans_assign(point.distances)[1L]
        origin <- labels[i]
        if (destination != origin && counts[origin] > 1L) {
          labels[i] <- destination
          counts[origin] <- counts[origin] - 1L
          counts[destination] <- counts[destination] + 1L
          centers[[origin]] <- riem_kmeans_mean(riemobj, which(labels == origin), geometry,
                                               mean.maxiter, mean.eps)
          centers[[destination]] <- riem_kmeans_mean(riemobj, which(labels == destination),
                                                    geometry, mean.maxiter, mean.eps)
        }
      }
    }
    distances <- riem_kmeans_distances(riemobj$data, centers, geometry)
    final.labels <- riem_kmeans_assign(distances)
    objective <- riem_kmeans_objective(distances, final.labels)
    if (algorithm == "lloyd" && objective > old.objective + eps * abs(old.objective)) {
      centers <- old.centers
      distances <- old.distances
      labels <- old.labels
      objective <- old.objective
      termination <- "objective_increase_rejected"
      break
    }
    history <- c(history, objective)
    if (identical(final.labels, labels) && !any(tabulate(final.labels, k) == 0L)) {
      labels <- final.labels
      converged <- TRUE
      termination <- "stable_assignment"
      break
    }
    labels <- final.labels
  }
  # Returned assignments and score always refer to the final stored centers.
  labels <- riem_kmeans_assign(distances)
  objective <- riem_kmeans_objective(distances, labels)
  list(cluster = labels, centers = centers, score = objective,
       converged = converged, termination = termination, iterations = iteration,
       objective_history = history, initialization = seeds, empty_repairs = empty.repairs,
       empty_clusters = which(tabulate(labels, k) == 0L))
}

riem_kmeans_check_fit <- function(object) {
  if (!inherits(object, "riem_kmeans") || !identical(object$schema_version, 1L) ||
      !is.list(object$centers) || is.null(object$geometry) || is.null(object$input_template)) {
    stop("This is not a supported fitted clustering object; refit with riem.kmeans().",
         call. = FALSE)
  }
  template <- object$input_template
  template$data <- object$centers
  riem_resolve_geometry(template, object$geometry, capability = "mean")
  invisible(TRUE)
}
