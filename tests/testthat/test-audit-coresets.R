coreset_reference_probability <- function(x) {
  delta <- sweep(x, 2, colMeans(x), "-")
  squared <- rowSums(delta^2)
  if (sum(squared) == 0) rep(1 / nrow(x), nrow(x)) else
    0.5 / nrow(x) + 0.5 * squared / sum(squared)
}

test_that("lightweight coresets use iid draws and matching importance weights", {
  x <- matrix(c(0, 1, 10), ncol = 1)
  X <- wrap.euclidean(x)
  q <- coreset_reference_probability(x)
  set.seed(71)
  reference <- sample.int(nrow(x), 12, replace = TRUE, prob = q)
  set.seed(71)
  result <- riem.coreset18B(X, M = 12)
  expect_identical(result$coreid, reference)
  expect_true(anyDuplicated(result$coreid) > 0)
  expect_equal(result$weight, 1 / (12 * q[reference]), tolerance = 1e-14)
  native <- Riemann:::learning_coreset18B("euclidean", "intrinsic", X$data, 2L, 50L, 1e-5)
  expect_equal(as.vector(native$qx), q, tolerance = 1e-14)

  # Enumerate the nine iid ordered draws, including repeated indices. This is
  # an exact importance-sampling identity, not a Monte Carlo tolerance test.
  outcomes <- as.matrix(expand.grid(first = 1:3, second = 1:3))
  for (centers in list(0, c(-2, 4), c(0, 1, 10))) {
    costs <- apply(outer(x[, 1], centers, "-")^2, 1, min)
    probabilities <- apply(outcomes, 1, function(ids) prod(q[ids]))
    estimates <- apply(outcomes, 1, function(ids) sum(costs[ids] / (2 * q[ids])))
    expect_equal(sum(probabilities * estimates), sum(costs), tolerance = 1e-12)
  }
  set.seed(71)
  expect_identical(riem.coreset18B(X, M = 12), result)
  expect_length(riem.coreset18B(X)$coreid, 2L)
})

test_that("zero-dispersion coresets have a defined uniform design", {
  X <- wrap.euclidean(matrix(rep(3, 4), ncol = 1))
  set.seed(11)
  result <- riem.coreset18B(X, M = 9)
  expect_equal(result$weight, rep(4 / 9, 9))
  expect_equal(sum(result$weight * 3^2), 4 * 3^2)
  fit <- riem.kmeans18B(X, k = 1, M = 9, nstart = 1)
  expect_equal(as.vector(fit$means), 3)
  expect_equal(fit$score, 0)
  expect_true(fit$converged)
  expect_warning(tied <- riem.kmeans18B(X, k = 2, M = 9, nstart = 1), "empty cluster")
  expect_equal(dim(tied$means), c(1L, 1L, 2L))
  expect_equal(tied$empty_clusters, 2L)
  expect_equal(tied$cluster, rep(1L, 4))
})

test_that("one-cluster coreset fitting solves the weighted Euclidean problem", {
  x <- matrix(c(0, 1, 10), ncol = 1)
  q <- coreset_reference_probability(x)
  X <- wrap.euclidean(x)
  set.seed(1)
  ids <- sample.int(3, 12, replace = TRUE, prob = q)
  weights <- 1 / (12 * q[ids])
  expected <- sum(weights * x[ids, 1]) / sum(weights)
  set.seed(1)
  fit <- riem.kmeans18B(X, k = 1, M = 12, nstart = 1)
  expect_identical(fit$coreset$coreid, ids)
  expect_equal(fit$coreset$weight, weights, tolerance = 1e-14)
  expect_equal(as.vector(fit$means), expected, tolerance = 1e-12)
  expect_equal(fit$score, sum((x[, 1] - expected)^2), tolerance = 1e-12)
  expect_equal(tail(fit$objective_history, 1), sum(weights * (x[ids, 1] - expected)^2),
               tolerance = 1e-12)
  expect_true(all(diff(fit$objective_history) <= 1e-10))
  expect_equal(fit$iterations, 1L)
  expect_true(fit$converged)
})

coreset_reference_first_lloyd_pass <- function(x, k, M) {
  q <- coreset_reference_probability(x)
  ids <- sample.int(nrow(x), M, replace = TRUE, prob = q)
  sampled <- x[ids, , drop = FALSE]
  weights <- 1 / (M * q[ids])
  chosen <- integer(k)
  chosen[1] <- sample.int(M, 1, prob = weights)
  distance <- function(centers) vapply(seq_len(nrow(centers)), function(j)
    sqrt(rowSums(sweep(sampled, 2, centers[j, ], "-")^2)), numeric(M))
  nearest <- distance(sampled[chosen[1], , drop = FALSE])[, 1]
  for (j in 2:k) {
    chosen[j] <- if (max(nearest) == 0) chosen[1] else
      sample.int(M, 1, prob = weights * (nearest / max(nearest))^2)
    nearest <- pmin(nearest, distance(sampled[chosen[j], , drop = FALSE])[, 1])
  }
  centers <- sampled[chosen, , drop = FALSE]
  labels <- max.col(-distance(centers), ties.method = "first")
  for (j in seq_len(k)) {
    inside <- which(labels == j)
    if (length(inside)) centers[j, ] <- colSums(sampled[inside, , drop = FALSE] * weights[inside]) /
      sum(weights[inside])
  }
  centers
}

test_that("maxiter is honored and a Lloyd pass uses importance-weighted means", {
  x <- matrix(c(-8, -7, -6, -2, -1, 0, 1, 2, 5, 6, 11, 12), ncol = 1)
  X <- wrap.euclidean(x)
  set.seed(17)
  expected <- coreset_reference_first_lloyd_pass(x, 3L, 30L)
  set.seed(17)
  fit <- suppressWarnings(riem.kmeans18B(X, k = 3, M = 30, nstart = 1, maxiter = 1))
  expect_equal(as.vector(fit$means), as.vector(expected), tolerance = 1e-12)
  expect_equal(fit$iterations, 1L)
  expect_length(fit$objective_history, 2L)
  distances <- abs(outer(x[, 1], expected[, 1], "-"))
  labels <- max.col(-distances, ties.method = "first")
  expect_equal(fit$cluster, labels)
  expect_equal(fit$score, sum(distances[cbind(seq_along(labels), labels)]^2), tolerance = 1e-12)
})

test_that("log-Euclidean coreset means use the requested geometry and weights", {
  logs <- c(-2, 0.1, 1, 3)
  X <- wrap.spd(lapply(logs, function(x) matrix(exp(x), 1, 1)))
  set.seed(37)
  fit <- riem.kmeans18B(X, k = 1, M = 10, geometry = "log_euclidean", nstart = 1)
  expected <- weighted.mean(logs[fit$coreset$coreid], fit$coreset$weight)
  expect_equal(log(as.vector(fit$means)), expected, tolerance = 1e-10)
  expect_equal(fit$score, sum((logs - expected)^2), tolerance = 1e-10)
})

test_that("coreset starts are selected by their full-data squared-distance cost", {
  x <- matrix(c(0, 1, 10), ncol = 1)
  X <- wrap.euclidean(x)
  set.seed(218)
  starts <- replicate(5, riem.kmeans18B(X, k = 1, M = 2, nstart = 1), simplify = FALSE)
  scores <- vapply(starts, function(fit) sum((x[, 1] - as.vector(fit$means))^2), numeric(1))
  expected <- starts[[which.min(scores)]]
  set.seed(218)
  actual <- riem.kmeans18B(X, k = 1, M = 2, nstart = 5)
  expect_equal(actual$starts$score, scores, tolerance = 1e-12)
  expect_equal(actual$score, min(scores), tolerance = 1e-12)
  expect_identical(actual$coreset, expected$coreset)
  expect_equal(actual$means, expected$means, tolerance = 1e-12)
})

test_that("default coreset sizes work for singleton data and k equal to N", {
  singleton <- wrap.euclidean(matrix(7, 1, 1))
  expect_equal(riem.coreset18B(singleton), list(coreid = 1L, weight = 1))
  fit <- riem.kmeans18B(singleton, k = 1, nstart = 1)
  expect_equal(as.vector(fit$means), 7)
  expect_equal(fit$score, 0)
  X <- wrap.euclidean(matrix(c(0, 2), ncol = 1))
  set.seed(9)
  fit <- suppressWarnings(riem.kmeans18B(X, k = 2, nstart = 1))
  expect_length(fit$coreset$coreid, 2L)
  expect_equal(dim(fit$means), c(1L, 1L, 2L))
  expect_true(all(is.finite(fit$means)))
  expect_true(all(fit$cluster %in% 1:2))
})

test_that("coreset controls reject malformed values instead of silently coercing", {
  X <- wrap.euclidean(matrix(c(0, 1, 10), ncol = 1))
  for (bad in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 2), 1i)) {
    expect_error(riem.coreset18B(X, M = bad), "M")
    expect_error(riem.kmeans18B(X, k = 1, M = bad), "M")
    expect_error(riem.coreset18B(X, maxiter = bad), "maxiter")
    expect_error(riem.kmeans18B(X, maxiter = bad), "maxiter")
    expect_error(riem.kmeans18B(X, nstart = bad), "nstart")
  }
  expect_error(riem.kmeans18B(X, k = 2, M = 1), "M")
  expect_error(riem.kmeans18B(X, k = 4), "k")
  expect_error(riem.coreset18B(X, eps = 0), "eps")
  expect_error(riem.coreset18B(X, unused = 1), "Extra arguments")
  expect_error(riem.kmeans18B(X, unused = 1), "Extra arguments")
  expect_error(riem.coreset18B(list(data = list(1))), "riemdata")
})
