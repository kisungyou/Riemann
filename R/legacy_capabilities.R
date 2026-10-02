# Contracts shared by legacy ordination and inferential entry points.
riem_legacy_positive <- function(x, name) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x <= 0) stop(name, " must be a finite positive number.", call. = FALSE)
  as.double(x)
}

riem_legacy_dimensions <- function(ndim, n) {
  if (n < 2L) stop("Ordination requires at least two observations.", call. = FALSE)
  riem_regression_integer(ndim, "ndim", 1L, n - 1L)
}

riem_legacy_distances <- function(riemobj, geometry = NULL) {
  spec <- riem_resolve_geometry(riemobj, geometry, capability = "distance")
  distances <- basic_pdist(riemobj$name, riemobj$data, spec$backend)
  if (!is.matrix(distances) || any(!is.finite(distances)) || any(distances < 0)) {
    stop("The selected geometry returned invalid distances.", call. = FALSE)
  }
  list(distances = distances, geometry = spec)
}

riem_legacy_eigen <- function(gram, ndim, negative = c("truncate", "error"), absolute_scale = 0) {
  negative <- match.arg(negative)
  fit <- eigen((gram + t(gram)) / 2, symmetric = TRUE)
  scale <- max(abs(fit$values), absolute_scale)
  threshold <- if (scale == 0) 0 else 100 * .Machine$double.eps * nrow(gram) * scale
  bad <- fit$values < -threshold
  if (any(bad) && negative == "error") {
    stop("A materially negative eigenvalue prevents a positive-semidefinite representation.", call. = FALSE)
  }
  positive <- which(fit$values > threshold)
  used <- utils::head(positive, ndim)
  embedding <- matrix(0, nrow(gram), ndim)
  if (length(used)) embedding[, seq_along(used)] <-
    sweep(fit$vectors[, used, drop = FALSE], 2L, sqrt(fit$values[used]), "*")
  list(embed = embedding, eigenvalues = fit$values, effective_dimension = length(used),
       rank = length(positive), negative_eigenvalues = fit$values[bad],
       negative_policy = negative, eigen_tolerance = threshold)
}

riem_legacy_cmds <- function(distances, ndim, negative = c("truncate", "error")) {
  negative <- match.arg(negative)
  scale <- max(distances)
  d <- if (scale == 0) distances else distances / scale
  squared <- d^2
  gram <- -0.5 * (sweep(sweep(squared, 1L, rowMeans(squared), "-"),
                        2L, colMeans(squared), "-") + mean(squared))
  out <- riem_legacy_eigen(gram, ndim, negative)
  out$embed <- out$embed * scale
  out$eigenvalues <- out$eigenvalues * scale * scale
  out$negative_eigenvalues <- out$negative_eigenvalues * scale * scale
  out$eigen_tolerance <- out$eigen_tolerance * scale * scale
  if (any(!is.finite(out$embed)) || any(!is.finite(out$eigenvalues))) {
    stop("The ordination spectrum is not representable; rescale the input units explicitly.", call. = FALSE)
  }
  reconstructed <- as.matrix(stats::dist(if (scale == 0) out$embed else out$embed / scale))
  upper <- upper.tri(d)
  denominator <- sum(d[upper]^2)
  out$stress <- if (denominator == 0) 0 else
    sqrt(sum((d[upper] - reconstructed[upper])^2) / denominator)
  out$requested_dimension <- ndim
  out
}

riem_mutual_knn_graph <- function(distances, nnbd) {
  n <- nrow(distances)
  adjacency <- matrix(FALSE, n, n)
  for (i in seq_len(n)) {
    candidates <- setdiff(seq_len(n), i)
    neighbors <- candidates[order(distances[i, candidates], candidates)][seq_len(nnbd)]
    adjacency[i, neighbors] <- TRUE
  }
  adjacency <- adjacency & t(adjacency)
  paths <- matrix(Inf, n, n)
  paths[adjacency] <- distances[adjacency]
  diag(paths) <- 0
  for (k in seq_len(n)) paths <- pmin(paths, outer(paths[, k], paths[k, ], "+"))
  component <- integer(n)
  for (i in seq_len(n)) if (component[i] == 0L) {
    component[is.finite(paths[i, ])] <- max(component) + 1L
  }
  list(distance = paths, adjacency = adjacency, component = component)
}

riem_legacy_parameters <- function(parameters, allowed) {
  if (length(parameters)) {
    supplied <- names(parameters)
    if (is.null(supplied) || anyNA(supplied) || any(!nzchar(supplied)) ||
        anyDuplicated(supplied) || any(!supplied %in% allowed)) {
      stop("Extra arguments must have unique names from: ", paste(allowed, collapse = ", "), call. = FALSE)
    }
  }
  parameters
}

riem_transport_inputs <- function(x, y, p, geometry) {
  riem_check_newdata(x, y)
  spec <- riem_resolve_geometry(x, geometry, "distance")
  other <- riem_resolve_geometry(y, geometry, "distance")
  if (!identical(spec, other)) stop("Both samples must use the same geometry.", call. = FALSE)
  p <- riem_legacy_positive(p, "p")
  if (p < 1) stop("p must be at least one for a Wasserstein distance.", call. = FALSE)
  list(geometry = spec, p = p)
}

riem_transport_solve <- function(distances, p, wx, wy, smooth = FALSE) {
  if (!is.matrix(distances) || any(!is.finite(distances)) || any(distances < 0)) {
    stop("Transport requires a finite nonnegative distance matrix.", call. = FALSE)
  }
  scale <- max(distances)
  normalized <- if (scale > 0) distances / scale else distances
  if (scale == 0 || nrow(distances) == 1L || ncol(distances) == 1L) {
    output <- list(plan = outer(wx, wy))
  } else {
    output <- if (smooth) T4transport_ipotD(normalized, p, wx, wy) else
      T4transport::wassersteinD(normalized, p, wx = wx, wy = wy)
  }
  plan <- output$plan
  if (!is.matrix(plan) || !identical(dim(plan), dim(distances)) ||
      any(!is.finite(plan)) || any(plan < -1e-12)) {
    stop("The transport solver returned an invalid plan.", call. = FALSE)
  }
  residual <- max(abs(rowSums(plan) - wx), abs(colSums(plan) - wy))
  if (residual > 1e-6) stop("The transport plan does not satisfy its marginal constraints.", call. = FALSE)
  plan <- pmax(plan, 0)
  output$plan <- plan
  output$distance <- scale * sum(plan * normalized^p)^(1 / p)
  if (!is.finite(output$distance)) stop("The transport distance is not numerically representable.", call. = FALSE)
  output$marginal_residual <- residual
  output$solver <- if (smooth) "IPOT_approximation" else "linear_program"
  output$validation_status <- if (smooth) "experimental_approximation_no_optimality_certificate" else
    "finite_plan_and_marginals_checked"
  output
}

riem_angular_parameter <- function(parameter, name) {
  parameter <- as.matrix(parameter)
  if (!check_spdmat(parameter)) stop(name, " must be a finite symmetric positive-definite matrix.", call. = FALSE)
  scaled <- parameter / max(abs(parameter))
  scaled <- scaled * (nrow(scaled) / sum(diag(scaled)))
  if (!check_spdmat(scaled)) stop(name, " is numerically singular after scale normalization.", call. = FALSE)
  scaled
}

riem_angular_density <- function(data, parameter, log) {
  if (!is.logical(log) || length(log) != 1L || is.na(log)) stop("log must be TRUE or FALSE.", call. = FALSE)
  p <- nrow(parameter)
  r <- ncol(data[[1L]])
  if (nrow(data[[1L]]) != p) stop("The parameter matrix must match the observations' row dimension.", call. = FALSE)
  root <- chol(parameter)
  logdet <- 2 * sum(base::log(diag(root)))
  values <- vapply(data, function(x) {
    whitened <- forwardsolve(t(root), x)
    smallroot <- chol(crossprod(whitened))
    -r * logdet / 2 - p * sum(base::log(diag(smallroot)))
  }, numeric(1))
  if (log) values else exp(values)
}

riem_angular_mle <- function(data, parameters, matrix_variant) {
  parameters <- riem_legacy_parameters(parameters, c("maxiter", "eps"))
  maxiter <- if (is.null(parameters$maxiter)) 50L else
    riem_regression_integer(parameters$maxiter, "maxiter", 1L)
  eps <- if (is.null(parameters$eps)) 1e-5 else riem_legacy_positive(parameters$eps, "eps")
  p <- nrow(data[[1L]])
  r <- ncol(data[[1L]])
  if (matrix_variant && r == p && p > 1L) {
    stop("When r equals p the MACG density is uniform for every Sigma; its shape parameter is not identifiable.", call. = FALSE)
  }
  span <- Reduce(`+`, lapply(data, tcrossprod))
  if (!check_spdmat(span)) stop("The observations do not span the ambient space; no positive-definite shape estimate is available.", call. = FALSE)
  estimate <- diag(p)
  converged <- FALSE
  for (iteration in seq_len(maxiter)) {
    root <- chol(estimate)
    update <- Reduce(`+`, lapply(data, function(x) {
      whitened <- forwardsolve(t(root), x)
      cross <- crossprod(whitened)
      x %*% solve(cross, t(x))
    }))
    update <- riem_angular_parameter(update, "Shape update")
    step <- norm(update - estimate, "F")
    estimate <- update
    if (step <= eps) {
      converged <- TRUE
      break
    }
  }
  attr(estimate, "converged") <- converged
  attr(estimate, "iterations") <- iteration
  attr(estimate, "termination") <- if (converged) "iterate_tolerance" else "maxiter"
  attr(estimate, "step_norm") <- step
  attr(estimate, "objective") <- -mean(riem_angular_density(data, estimate, TRUE))
  attr(estimate, "validation_status") <- "fixed_point_stationarity_only_existence_and_uniqueness_not_certified"
  if (!converged) warning("The angular shape iteration did not converge; the last positive-definite estimate is returned.", call. = FALSE)
  estimate
}
