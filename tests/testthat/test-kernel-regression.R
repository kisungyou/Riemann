# Independent references use only Euclidean distances or the sphere chord.
reference_kernel <- function(distance, response, bandwidth) {
  vapply(seq_len(ncol(distance)), function(j) {
    logweight <- -(distance[, j] / bandwidth)^2 / 2
    weight <- exp(logweight - max(logweight))
    sum(weight * response) / sum(weight)
  }, numeric(1))
}

reference_cv <- function(x, y, bandwidth, foldid) {
  out <- numeric(length(y))
  for (fold in unique(foldid)) {
    test <- which(foldid == fold)
    train <- which(foldid != fold)
    distance <- abs(outer(x[train], x[test], "-"))
    out[test] <- reference_kernel(distance, y[train], bandwidth)
  }
  out
}

test_that("the Gaussian smoother matches an independent reference", {
  x <- c(-2, -0.3, 0.2, 1.7, 4)
  y <- c(2, -1, 3, 0, 5)
  X <- wrap.euclidean(matrix(x, ncol = 1))
  fit <- riem.m2skreg(X, y, bandwidth = 0.8)
  expected <- reference_kernel(abs(outer(x, x, "-")), y, 0.8)
  expect_equal(fit$ypred, expected, tolerance = 1e-12)
  expect_equal(fitted(fit), expected, tolerance = 1e-12)
  expect_equal(residuals(fit), y - expected, tolerance = 1e-12)
  expect_equal(summary(fit)$training_rmse, sqrt(mean((y - expected)^2)))
  expect_identical(fit$geometry$geometry_id, "euclidean")
  expect_equal(predict(fit, X), expected, tolerance = 1e-12)
  expect_output(print(fit), "Bandwidth")
  expect_output(print(summary(fit)), "self-weights")
})

test_that("extrinsic sphere geometry survives fitting and prediction", {
  theta <- c(0, 0.4, 1.2, 2.2)
  x <- cbind(cos(theta), sin(theta))
  y <- c(0, 1, -2, 3)
  X <- wrap.sphere(x)
  fit <- riem.m2skreg(X, y, bandwidth = 0.8, geometry = "extrinsic")
  expected <- reference_kernel(as.matrix(dist(x)), y, 0.8)
  expect_identical(fit$geometry$geometry_id, "chordal")
  expect_equal(predict(fit, X), expected, tolerance = 1e-10)
  expect_equal(predict(fit, X, geometry = "chordal"), expected, tolerance = 1e-10)
  expect_error(predict(fit, X, geometry = "intrinsic"), "must match")
  intrinsic <- riem.m2skreg(X, y, bandwidth = 0.8, geometry = "intrinsic")
  expect_gt(max(abs(intrinsic$ypred - fit$ypred)), 1e-3)
})

test_that("prediction is invariant to batching, blocking and serialization", {
  X <- wrap.euclidean(matrix(c(0, 1, 3), ncol = 1))
  fit <- riem.m2skreg(X, c(1, -2, 4), bandwidth = 0.4)
  grid <- wrap.euclidean(matrix(c(0.5, 2, 10), ncol = 1))
  expected <- predict(fit, grid)
  expect_equal(predict(fit, grid, block_size = 1), expected)
  expect_equal(predict(fit, wrap.euclidean(list(2))), expected[2])
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_equal(predict(readRDS(path), grid), expected)
  supported <- predict(fit, grid, diagnostics = TRUE)
  expect_equal(supported$prediction, expected)
  expect_equal(supported$diagnostics$nearest_distance, c(0.5, 1, 7))
  expect_true(all(supported$diagnostics$effective_n >= 1))
  expect_true(all(supported$diagnostics$effective_n <= 3 + 1e-12))
})

test_that("legacy objects require an explicit geometry", {
  X <- wrap.euclidean(matrix(0:2, ncol = 1))
  fit <- riem.m2skreg(X, 0:2)
  legacy <- fit
  legacy$geometry <- NULL
  legacy$schema_version <- NULL
  expect_error(predict(legacy, X), "legacy fit has no geometry")
  expect_warning(value <- predict(legacy, X, geometry = "intrinsic"), "cannot be verified")
  expect_equal(value, fit$ypred)
  expect_null(legacy$geometry)
  bad_schema <- fit
  bad_schema$schema_version <- 99L
  expect_error(predict(bad_schema, X), "Unsupported.*schema")
})

test_that("invalid response, bandwidth and prediction arguments are rejected", {
  X <- wrap.euclidean(matrix(0:2, ncol = 1))
  for (h in list(0, -1, Inf, NA_real_, numeric(), c(1, 2), "1", 1 + 1i)) {
    expect_error(riem.m2skreg(X, 0:2, bandwidth = h), "bandwidth")
  }
  for (y in list(c(1, 2), c(1, NA, 3), c(1, Inf, 3), letters[1:3],
                 matrix(1:3, ncol = 1), factor(1:3), c(1i, 2i, 3i))) {
    expect_error(riem.m2skreg(X, y), "response")
  }
  fit <- riem.m2skreg(X, 0:2)
  expect_error(predict(fit, X, block_size = 0), "block_size")
  expect_error(predict(fit, X, diagnostics = NA), "diagnostics")
  expect_error(predict(fit, X, unknown = TRUE), "Unknown")
  expect_error(predict(fit, wrap.euclidean(matrix(1:4, ncol = 2))), "incompatible")
  ordered <- X
  ordered$feature_names <- "training_feature"
  fit_ordered <- riem.m2skreg(ordered, 0:2)
  expect_error(predict(fit_ordered, X), "feature")
})

test_that("tiny bandwidth and huge distances retain finite relative weights", {
  weights <- Riemann:::riem_kernel_weights(c(1e300, 1e300, 2e300), 1e-300)
  expect_equal(weights, c(0.5, 0.5, 0))
  d <- matrix(c(1e300, 2e300, 3e300), ncol = 1)
  result <- Riemann:::riem_kernel_predict_distances(d, c(1, 3, 5), 1e300)
  expect_equal(result$prediction, reference_kernel(d / 1e300, c(1, 3, 5), 1),
               tolerance = 1e-12)
  tied <- Riemann:::riem_kernel_predict_distances(matrix(1e300, 2, 1),
                                                 c(-1e308, 1e308), 1e-300)
  expect_equal(tied$prediction, 0)
  expect_equal(tied$diagnostics$effective_n, 2)
  expect_true(is.infinite(tied$diagnostics$nearest_over_bandwidth))
  expect_error(Riemann:::riem_kernel_weights(c(0, Inf), 1), "finite")
  expect_error(Riemann:::riem_kernel_weights(c(-1, 0), 1), "nonnegative")
  X <- wrap.euclidean(matrix(c(0, 1, 3), ncol = 1))
  fit <- riem.m2skreg(X, c(0, 2, 4), bandwidth = 1e-100)
  expect_equal(predict(fit, wrap.euclidean(list(100))), 4)
})

test_that("CV aggregates every fold and preserves the complete error table", {
  x <- 0:5
  y <- c(0, 3, 0, 0, 0, 0)
  folds <- rep(1:3, each = 2)
  bandwidths <- c(0.1, 0.5, 1, 2, 10)
  X <- wrap.euclidean(matrix(x, ncol = 1))
  expected <- vapply(bandwidths, function(h) {
    sum((reference_cv(x, y, h, folds) - y)^2)
  }, numeric(1))
  fit <- riem.m2skregCV(X, y, bandwidths, foldid = folds)
  expect_equal(fit$errors[, "SSE"], expected, tolerance = 1e-12)
  expect_equal(dim(fit$errors), c(5L, 2L))
  expect_equal(dim(fit$fold_errors), c(5L, 3L))
  expect_equal(rowSums(fit$fold_errors), expected, tolerance = 1e-12)
  expect_equal(fit$bandwidth, 2)
  expect_equal(fit$cv_prediction, reference_cv(x, y, 2, folds), tolerance = 1e-12)
  expect_equal(fit$ypred, riem.m2skreg(X, y, 2)$ypred)
  expect_equal(predict(fit, X), fit$ypred)
  expect_identical(fit$foldid, folds)
  expect_true(all(fit$candidate_status == "ok"))
  expect_true(all(is.na(fit$fold_failure)))
})

test_that("CV supports single held-out observations and one training observation", {
  X <- wrap.euclidean(matrix(c(0, 2), ncol = 1))
  fit <- riem.m2skregCV(X, c(1, 4), bandwidths = c(1e-100, 1), kfold = 2)
  expect_equal(fit$errors[, "SSE"], c(18, 18))
  expect_equal(fit$cv_prediction, c(4, 1))
  expect_equal(fit$bandwidth, 1e-100)
  expect_equal(dim(fit$errors), c(2L, 2L))
  singleton <- riem.m2skreg(wrap.euclidean(list(2)), 7)
  expect_equal(singleton$ypred, 7)
  expect_equal(predict(singleton, wrap.euclidean(list(10))), 7)
  expect_error(riem.m2skregCV(wrap.euclidean(list(2)), 7), "at least two")
})

test_that("fold inputs are validated and generated folds are reproducible", {
  X <- wrap.euclidean(matrix(0:5, ncol = 1))
  for (k in list(0, 1, 2.2, 7, Inf, NA_real_)) {
    expect_error(riem.m2skregCV(X, 0:5, kfold = k), "kfold")
  }
  for (f in list(rep(1, 6), 1:3, c(1, 2, 3, 1, 2, NA),
                 c(1, 2, 3, 1, 2, Inf), c("a", "b", "a", "b", "a", ""))) {
    expect_error(riem.m2skregCV(X, 0:5, foldid = f), "foldid")
  }
  expect_error(riem.m2skregCV(X, 0:5, kfold = 2, foldid = rep(1:3, each = 2)),
               "does not match")
  expect_error(riem.m2skregCV(X, 0:5, bandwidths = numeric()), "bandwidths")
  set.seed(428)
  a <- riem.m2skregCV(X, 0:5, bandwidths = 1, kfold = 3)
  set.seed(428)
  b <- riem.m2skregCV(X, 0:5, bandwidths = 1, kfold = 3)
  expect_identical(a$foldid, b$foldid)
  expect_equal(a$errors, b$errors)
  expect_equal(as.vector(table(a$foldid)), rep(2L, 3))
})

test_that("CV records failed candidates and refuses all-failed tuning", {
  X <- wrap.euclidean(matrix(c(0, 0, 10, 10), ncol = 1))
  y <- c(0, 0, 1e200, 1e200)
  folds <- c("a", "b", "a", "b")
  fit <- riem.m2skregCV(X, y, bandwidths = c(0.01, 10), foldid = folds)
  expect_equal(fit$bandwidth, 0.01)
  expect_identical(fit$candidate_status, c("ok", "failed"))
  expect_equal(unname(fit$errors[1, "SSE"]), 0)
  expect_true(is.infinite(fit$errors[2, "SSE"]))
  expect_true(all(grepl("nonfinite", fit$fold_failure[2, ])))
  expect_equal(dim(fit$errors), c(2L, 2L))
  expect_error(riem.m2skregCV(X, y, bandwidths = 10, foldid = folds),
               "No bandwidth has finite errors")
})
