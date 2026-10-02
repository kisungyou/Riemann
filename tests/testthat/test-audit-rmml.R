rmml_reference_distances <- function(embedded, labels, lambda) {
  p <- ncol(embedded)
  same <- different <- matrix(0, p, p)
  for (i in seq_len(nrow(embedded) - 1L)) for (j in (i + 1L):nrow(embedded)) {
    scatter <- tcrossprod(embedded[i, ] - embedded[j, ])
    if (labels[i] == labels[j]) same <- same + scatter else different <- different + scatter
  }
  power <- function(x, exponent) {
    eig <- eigen(x, symmetric = TRUE)
    eig$vectors %*% diag(eig$values^exponent, nrow(x)) %*% t(eig$vectors)
  }
  lhs <- solve(same + lambda * diag(p))
  root <- power(lhs, 0.5)
  inverse <- power(lhs, -0.5)
  metric <- root %*% power(inverse %*% (different + lambda * diag(p)) %*% inverse, 0.5) %*% root
  as.matrix(dist(embedded %*% power(metric, 0.5)))
}

test_that("RMML handles projector embeddings of proper Grassmann subspaces", {
  points <- lapply(c(0, 0.2, 0.8, 1), function(t) matrix(c(cos(t), sin(t), 0), 3, 1))
  labels <- c(1, 1, 2, 2)
  embedded <- t(vapply(points, function(x) as.vector(tcrossprod(x)), numeric(9)))
  expected <- rmml_reference_distances(embedded, labels, 0.1)
  actual <- riem.rmml(wrap.grassmann(points), labels)
  expect_equal(unname(actual), unname(expected), tolerance = 1e-8)
  expect_true(all(is.finite(actual)))
  expect_equal(actual, t(actual), tolerance = 1e-14)
  expect_equal(diag(actual), numeric(4))
  reversed <- Map(function(x, sign) x * sign, points, c(-1, 1, -1, 1))
  expect_equal(riem.rmml(wrap.grassmann(reversed), labels), actual, tolerance = 1e-8)
  expect_s3_class(riem.rmml(wrap.grassmann(points), labels, as.dist = TRUE), "dist")
})

test_that("RMML matches an independent embedded Euclidean calculation", {
  points <- rbind(c(0, 0), c(1, 0.2), c(3, 1), c(4, 2))
  labels <- c("a", "a", "b", "b")
  expected <- rmml_reference_distances(points, labels, 0.2)
  actual <- riem.rmml(wrap.euclidean(points), labels, lambda = 0.2)
  expect_equal(unname(actual), unname(expected), tolerance = 1e-9)
  labels[2] <- NA_character_
  subset <- riem.rmml(wrap.euclidean(points), labels, lambda = 0.2)
  expect_equal(unname(subset), unname(rmml_reference_distances(points[-2, ], labels[-2], 0.2)),
               tolerance = 1e-9)
})

test_that("RMML is invariant under nontrivial Grassmann basis rotations", {
  angles <- c(0, 0.2, 0.7, 1.1)
  points <- lapply(angles, function(theta) cbind(
    c(cos(theta), 0, sin(theta), 0),
    c(0, cos(theta / 2), 0, sin(theta / 2))))
  labels <- c(1, 1, 2, 2)
  embedded <- t(vapply(points, function(x) as.vector(tcrossprod(x)), numeric(16)))
  expected <- rmml_reference_distances(embedded, labels, 0.3)
  actual <- riem.rmml(wrap.grassmann(points), labels, lambda = 0.3)
  expect_equal(unname(actual), unname(expected), tolerance = 1e-8)
  rotated <- Map(function(frame, theta) frame %*%
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2),
    points, c(0.7, -0.3, 0.8, 1.3))
  expect_equal(riem.rmml(wrap.grassmann(rotated), labels, lambda = 0.3), actual,
               tolerance = 1e-8)
})

test_that("RMML validates label and regularization controls before native dispatch", {
  X <- wrap.euclidean(matrix(1:4, ncol = 1))
  expect_error(riem.rmml(X, c(1, 2)), "one label per observation")
  expect_error(riem.rmml(X, list(1, 1, 2, 2)), "label")
  expect_error(riem.rmml(X, rep(NA_real_, 4)), "two nonmissing")
  expect_error(riem.rmml(X, c(1, NA, NA, NA)), "two nonmissing")
  for (bad in list(NA_real_, Inf, numeric(), c(0.1, 0.2), 1i))
    expect_error(riem.rmml(X, c(1, 1, 2, 2), lambda = bad), "lambda")
  expect_error(riem.rmml(X, c(1, 1, 2, 2), as.dist = NA), "as.dist")
  expect_error(riem.rmml(X, c(1, 1, 2, 2), as.dist = 1), "as.dist")
})
