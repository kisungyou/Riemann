audit_sammon_loss <- function(original, coordinates) {
  target <- as.numeric(stats::dist(original))
  fitted <- as.numeric(stats::dist(coordinates))
  sum((target - fitted)^2 / target) / sum(target)
}

test_that("Sammon mapping descends on the previously divergent small-scale fixture", {
  set.seed(17)
  x <- matrix(stats::rnorm(30), 10, 3) * .01
  initial <- stats::cmdscale(stats::dist(x), k = 2)
  iterations <- c(1L, 2L, 5L, 20L, 100L)
  fits <- lapply(iterations, function(n) {
    riem.sammon(wrap.euclidean(x), ndim = 2, maxiter = n, eps = 0)
  })
  losses <- c(audit_sammon_loss(x, initial), vapply(fits, function(fit) {
    expect_true(all(is.finite(fit$embed)))
    target <- as.numeric(stats::dist(x))
    fitted <- as.numeric(stats::dist(fit$embed))
    expect_equal(fit$stress, sqrt(sum((target - fitted)^2) / sum(target^2)),
                 tolerance = 1e-12)
    audit_sammon_loss(x, fit$embed)
  }, 0))
  expect_true(all(diff(losses) <= 1e-12))
  expect_lt(tail(losses, 1), losses[1] / 2)
})

test_that("the first Sammon update agrees with derivatives of its stated scalar loss", {
  set.seed(17)
  x <- matrix(stats::rnorm(30), 10, 3)
  distance.scale <- max(stats::dist(x))
  normalized <- x / distance.scale
  y <- stats::cmdscale(stats::dist(normalized), k = 2)
  objective <- function(z) audit_sammon_loss(normalized, z)
  baseline <- objective(y)
  # Obtain both derivatives from finite differences of the scalar objective;
  # this oracle does not use the production gradient or Hessian expressions.
  gradient <- hessian <- matrix(0, nrow(y), ncol(y))
  h <- 1e-5
  for (index in seq_along(y)) {
    plus <- minus <- y
    plus[index] <- plus[index] + h
    minus[index] <- minus[index] - h
    fplus <- objective(plus)
    fminus <- objective(minus)
    gradient[index] <- (fplus - fminus) / (2 * h)
    hessian[index] <- (fplus - 2 * baseline + fminus) / h^2
  }
  direction <- -gradient / pmax(abs(hessian),
    max(1e-12 * max(abs(hessian)), .Machine$double.eps))
  direction <- direction / max(1, max(abs(direction)))
  slope <- sum(gradient * direction)
  step <- .3
  for (attempt in seq_len(60)) {
    expected <- scale(y + step * direction, center = TRUE, scale = FALSE)
    if (objective(expected) < baseline &&
        objective(expected) <= baseline + 1e-4 * step * slope) break
    step <- step / 2
  }
  expect_lt(objective(expected), baseline)
  actual <- riem.sammon(wrap.euclidean(x), ndim = 2, maxiter = 1, eps = 0)
  # Eigenvector signs and order need not agree between the two eigensolvers.
  expect_equal(as.numeric(stats::dist(actual$embed)) / distance.scale,
               as.numeric(stats::dist(expected)), tolerance = 2e-6)
})

test_that("Sammon fitting is equivariant to a common change of units", {
  set.seed(17)
  x <- matrix(stats::rnorm(30), 10, 3)
  baseline <- riem.sammon(wrap.euclidean(x), maxiter = 100, eps = 1e-8)
  for (multiplier in c(1e-6, .01, 1e6)) {
    fit <- riem.sammon(wrap.euclidean(multiplier * x), maxiter = 100, eps = 1e-8)
    expect_equal(as.numeric(stats::dist(fit$embed)) / multiplier,
                 as.numeric(stats::dist(baseline$embed)), tolerance = 1e-9)
    expect_equal(fit$stress, baseline$stress, tolerance = 1e-10)
  }
})

test_that("Sammon handles projected collisions and non-Euclidean initialization", {
  # The minor-axis observations coincide in the initial one-dimensional MDS
  # representation, despite having positive original distance.
  x <- rbind(c(-2, 0), c(2, 0), c(0, .1), c(0, -.1))
  fit <- riem.sammon(wrap.euclidean(x), ndim = 1, maxiter = 100)
  expect_true(all(is.finite(fit$embed)))
  expect_true(all(stats::dist(fit$embed) > 0))
  expect_true(is.finite(fit$stress))
  sphere <- wrap.sphere(rbind(c(1, 0), c(0, 1), c(-1, 0), c(0, -1)))
  curved <- riem.sammon(sphere, ndim = 3)
  expect_equal(dim(curved$embed), c(4L, 3L))
  expect_true(all(is.finite(curved$embed)))
  expect_true(is.finite(curved$stress))
})

test_that("Sammon rejects zero original distances and malformed controls", {
  duplicate <- wrap.euclidean(matrix(c(0, 0, 1), ncol = 1))
  expect_error(riem.sammon(duplicate, ndim = 1), "strictly positive.*coincident")
  sphere <- wrap.sphere(rbind(c(1, 0), c(2, 0), c(0, 1)))
  expect_error(riem.sammon(sphere, ndim = 1), "strictly positive.*coincident")
  x <- wrap.euclidean(matrix(c(0, 1, 3), ncol = 1))
  for (bad in list(0, -1, 1.5, Inf, NA, c(1, 2))) {
    expect_error(riem.sammon(x, ndim = bad), "ndim")
    expect_error(riem.sammon(x, maxiter = bad), "maxiter")
  }
  expect_error(riem.sammon(x, ndim = 3), "smaller")
  for (bad in list(-1, Inf, NA, c(1, 2), 1i)) {
    expect_error(riem.sammon(x, eps = bad), "eps")
  }
  expect_error(riem.sammon(x, typo = 1), "Extra arguments")
})

test_that("tangent PCA retains the analytic two-point variance for extreme weights", {
  x <- wrap.euclidean(matrix(c(0, 1), ncol = 1))
  # For any two strictly positive probability weights, the unbiased weighted
  # variance is (x2-x1)^2/2, irrespective of their imbalance.
  for (small in c(1e-12, 1e-16, 1e-100, 1e-300, 1e-320)) {
    for (weights in list(c(1, small), c(small, 1))) {
      fit <- riem.pga(x, ndim = 1, weight = weights)
      expect_equal(fit$variance, .5, tolerance = 1e-10)
      expect_equal(fit$total.variance, .5, tolerance = 1e-10)
      expect_equal(fit$rank, 1)
      expect_equal(fit$explained.variance, 1, tolerance = 1e-12)
    }
  }
})

test_that("tangent PCA weight rescaling and degenerate cases are explicit", {
  x <- wrap.euclidean(matrix(c(0, 1), ncol = 1))
  fit <- riem.pga(x, ndim = 1, weight = c(1e300, 1e284))
  expect_equal(fit$variance, .5, tolerance = 1e-12)
  expect_error(riem.pga(x, ndim = 1, weight = c(1e300, 1e-300)),
               "Positive weights underflowed")
  expect_warning(single <- riem.pga(x, ndim = 1, weight = c(1, 0)),
                 "numerical rank is 0")
  expect_equal(single$rank, 0)
  expect_equal(single$total.variance, 0)
  expect_equal(single$variance, numeric(0))
  expect_error(riem.pga(x, ndim = 1, weight = c(1, 0), center.tangent = FALSE),
               "at least two positive")
})

test_that("stationary mean diagnostics do not imply a local minimum", {
  angles <- c(.1, -.1, 3, -3)
  points <- cbind(sin(angles), 0, cos(angles))
  weights <- c(.3, .3, .2, .2)
  fit <- riem.mean(wrap.sphere(points), weight = weights, eps = 1e-12)
  expect_identical(fit$termination, "stationary")
  expect_true(fit$converged)
  expect_lt(fit$gradient_norm, 1e-12)
  expect_equal(as.numeric(fit$mean), c(0, 0, 1), tolerance = 1e-12)
  # An independent spherical-distance calculation finds a strict decrease in
  # a perpendicular tangent direction. The documented certificate is only
  # first-order stationarity, even when convergence is TRUE.
  perturbed <- c(0, sin(.001), cos(.001))
  smaller <- sum(weights * acos(pmax(-1, pmin(1, drop(points %*% perturbed))))^2)
  expect_lt(smaller, fit$objective - 1e-6)
  expect_lt(2 * sum(weights * angles / tan(angles)), 0)
})
