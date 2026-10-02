test_that("Lloyd clustering matches Euclidean WCSS and final-center assignment", {
  z <- rbind(c(-3, 0), c(-2, 0.3), c(-2.5, -0.2),
             c(3, 0), c(2, 0.2), c(2.5, -0.3))
  x <- wrap.euclidean(z)
  fit <- riem.kmeans(x, k = 2, init = c(1, 4), nstart = 1)
  expect_true(fit$converged)
  expected <- stats::kmeans(z, centers = z[c(1, 4), ], algorithm = "Lloyd")
  expect_equal(fit$score, expected$tot.withinss, tolerance = 1e-9)
  centers <- do.call(rbind, lapply(fit$centers, as.numeric))
  d2 <- vapply(seq_len(nrow(centers)), function(j) {
    rowSums(sweep(z, 2, centers[j, ])^2)
  }, numeric(nrow(z)))
  labels <- max.col(-d2, ties.method = "first")
  expect_identical(fit$cluster, labels)
  expect_equal(fit$score, sum(d2[cbind(seq_len(nrow(z)), labels)]), tolerance = 1e-10)
  expect_equal(predict(fit, x), fit$cluster)
  expect_equal(dim(fit$means), c(2, 1, 2))
  expect_true(all(diff(fit$objective_history) <= 1e-10))
  expect_null(fit$input_template$data)
})

test_that("k-means++ samples squared rather than raw distances", {
  x <- wrap.euclidean(matrix(c(0, 1, 4, 10), ncol = 1))
  geom <- riem.geometry(x)
  any.different <- FALSE
  for (seed in 1:30) {
    set.seed(seed)
    first <- sample.int(4, 1)
    distances <- abs(c(0, 1, 4, 10) - c(0, 1, 4, 10)[first])
    second <- sample.int(4, 1, prob = distances^2)
    set.seed(seed)
    actual <- riem_kmeans_initialize(x$data, 2, "plus", geom)
    expect_identical(actual, c(first, second))
    set.seed(seed)
    sample.int(4, 1)
    raw.second <- sample.int(4, 1, prob = distances)
    any.different <- any.different || raw.second != second
  }
  expect_true(any.different)
})

test_that("multiple starts and MacQueen are reproducible and objective-consistent", {
  x <- wrap.euclidean(rbind(c(-2, 0), c(-1, 1), c(0, 0), c(2, 1), c(3, 0)))
  set.seed(47)
  fit <- riem.kmeans(x, k = 2, nstart = 3, algorithm = "MacQueen")
  set.seed(47)
  repeat.fit <- riem.kmeans(x, k = 2, nstart = 3, algorithm = "MacQueen")
  expect_identical(fit$cluster, repeat.fit$cluster)
  expect_equal(fit$starts, repeat.fit$starts)
  expect_equal(nrow(fit$starts), 3)
  expect_equal(fit$score, min(fit$starts$objective[fit$starts$valid]))
  expect_identical(predict(fit, x), fit$cluster)
  d <- predict(fit, x, type = "distance")
  expect_equal(fit$score, sum(d[cbind(seq_along(fit$cluster), fit$cluster)]^2),
               tolerance = 1e-10)
})

test_that("clustering retains both SPD geometries through serialization and singleton prediction", {
  x <- wrap.spd(lapply(c(1, 1.5, 4, 6), function(a) diag(c(a, 1 / a))))
  for (geometry in c("affine_invariant", "log_euclidean")) {
    fit <- riem.kmeans(x, k = 2, geometry = geometry, init = c(1, 4), nstart = 1,
                       mean.eps = 1e-10)
    expect_true(fit$converged)
    expect_identical(fit$geometry$geometry_id, geometry)
    distances <- predict(fit, x, type = "distance")
    independent <- vapply(fit$centers, function(center) {
      vapply(x$data, function(a) sqrt(sum(log(diag(a) / diag(center))^2)), numeric(1))
    }, numeric(length(x$data)))
    expect_equal(unname(distances), independent, tolerance = 1e-8)
    singleton <- vapply(x$data, function(a) predict(fit, wrap.spd(list(a))), integer(1))
    expect_identical(singleton, fit$cluster)
    path <- tempfile(fileext = ".rds")
    saveRDS(fit, path)
    restored <- readRDS(path)
    unlink(path)
    expect_equal(predict(restored, x), fit$cluster)
  }
})

test_that("empty clusters are repaired and impossible ties fail visibly", {
  x <- wrap.euclidean(matrix(c(0, 0, 5, 5, 10, 10), ncol = 1))
  fit <- riem.kmeans(x, k = 3, init = c(1, 2, 3), nstart = 1)
  expect_true(fit$converged)
  expect_true(fit$diagnostics$empty_repairs > 0)
  expect_equal(fit$score, 0, tolerance = 1e-12)
  expect_equal(tabulate(fit$cluster, 3), rep(2L, 3))
  impossible <- wrap.euclidean(matrix(rep(1, 4), ncol = 1))
  error <- tryCatch(riem.kmeans(impossible, k = 2, nstart = 2), error = identity)
  expect_s3_class(error, "riem_kmeans_error")
  expect_true(all(!error$starts$valid))
  expect_match(conditionMessage(error), "distinct")
  tied <- riem.kmeans(wrap.euclidean(matrix(c(-1, 1), ncol = 1)),
                      k = 2, init = c(1, 2), nstart = 1)
  expect_identical(predict(tied, wrap.euclidean(matrix(0, 1, 1))), 1L)
})

test_that("clustering validates controls instead of silently rounding or clamping", {
  x <- wrap.euclidean(matrix(1:8, 4, 2))
  expect_error(riem.kmeans(x, k = 1.5), "integer")
  expect_error(riem.kmeans(x, k = 5), "integer")
  expect_error(riem.kmeans(x, maxiter = 0), "integer")
  expect_error(riem.kmeans(x, nstart = 0), "integer")
  expect_error(riem.kmeans(x, eps = 0), "positive")
  expect_error(riem.kmeans(x, init = c(1, 1)), "distinct")
  expect_error(riem.kmeans(x, nonexistent = 1), "control")
  fit <- riem.kmeans(x, k = 1, maxiter = 1, nstart = 1)
  expect_equal(fit$controls$maxiter, 1)
  expect_equal(fit$iterations, 1)
  expect_true(fit$converged)
  expect_error(predict(fit, wrap.sphere(matrix(c(1, 0), 1))), "incompatible")
})
