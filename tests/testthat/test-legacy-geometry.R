test_that("Grassmann principal angles retain small separations and quotient invariance", {
  X <- diag(4)[, 1:2]
  angles <- c(.2, .35)
  Y <- rbind(diag(cos(angles)), diag(sin(angles)))
  zero <- 0 * X
  first <- geometry_operations("grassmann", X, Y, zero, zero)
  maps <- geometry_operations("grassmann", X, Y, first$log, first$log)
  expect_equal(first$distance, sqrt(sum(angles^2)), tolerance = 1e-12)
  expect_equal(maps$metric, first$distance^2, tolerance = 1e-12)
  expect_equal(tcrossprod(maps$exp), tcrossprod(Y), tolerance = 1e-11)
  expect_equal(crossprod(X, first$log), matrix(0, 2, 2), tolerance = 1e-12)
  Q <- matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2)
  P <- diag(c(-1, 1))
  changed <- geometry_operations("grassmann", X %*% Q, Y %*% P, zero, zero)
  expect_equal(changed$distance, first$distance, tolerance = 1e-12)
  expect_equal(changed$log, first$log %*% Q, tolerance = 1e-11)

  tiny <- wrap.grassmann(list(matrix(c(1, 0, 0), 3),
    matrix(c(cos(1e-9), sin(1e-9), 0), 3)))
  expect_equal(riem.pdist(tiny)[1, 2] / 1e-9, 1, tolerance = 1e-6)
})

test_that("Grassmann cut loci and inverse-projection ambiguity are explicit", {
  X <- matrix(c(1, 0, 0), 3)
  Y <- matrix(c(0, 1, 0), 3)
  x <- wrap.grassmann(list(X, Y))
  expect_equal(riem.pdist(x)[1, 2], pi / 2)
  expect_error(geometry_operations("grassmann", X, Y, 0 * X, 0 * X), "nonunique")
  expect_error(riem.mean(x, geometry = "extrinsic"), "nonunique")
})

test_that("unsupported Stiefel and correlation intrinsic routes cannot bypass dispatch", {
  X <- diag(3)[, 1:2]
  stiefel <- wrap.stiefel(list(X, X))
  expect_error(riem.pdist(stiefel), "support|unavailable|restricted")
  expect_error(riem.mean(stiefel), "support|unavailable|restricted")
  expect_error(geometry_operations("stiefel", X, X, 0 * X, 0 * X), "unavailable")
  expect_equal(riem.pdist(stiefel, geometry = "extrinsic"), matrix(0, 2, 2))
  expect_error(riem.mean(wrap.stiefel(list(X, -X)), geometry = "extrinsic"), "nonunique")

  corr <- wrap.correlation(list(diag(2), matrix(c(1, .2, .2, 1), 2)))
  expect_error(riem.pdist(corr), "support|unavailable|restricted")
  expect_error(riem.mean(corr), "support|unavailable|restricted")
  expect_error(geometry_operations("correlation", diag(2), diag(2),
    matrix(0, 2, 2), matrix(0, 2, 2)), "unavailable")
})

test_that("rotation distances and body-coordinate logarithms have the declared scale", {
  rotate <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  X <- rotate(.2)
  Y <- rotate(.7)
  first <- geometry_operations("rotation", X, Y, 0 * X, 0 * X)
  maps <- geometry_operations("rotation", X, Y, first$log, first$log)
  expect_equal(first$distance, sqrt(2) * .5, tolerance = 1e-12)
  expect_equal(first$log + t(first$log), 0 * X, tolerance = 1e-12)
  expect_equal(maps$exp, Y, tolerance = 1e-12)
  expect_equal(maps$metric, first$distance^2, tolerance = 1e-12)
  expect_equal(riem.pdist(wrap.rotation(list(diag(2), rotate(pi))))[1, 2],
                sqrt(2) * pi, tolerance = 1e-12)
  expect_error(geometry_operations("rotation", diag(2), rotate(pi), 0 * X, 0 * X), "nonunique")
  expect_equal(riem.pdist(wrap.rotation(list(diag(2), rotate(1e-9))))[1, 2] /
    (sqrt(2) * 1e-9), 1, tolerance = 1e-6)
})

test_that("rotation projections stay in SO and common-arc intrinsic means agree", {
  data <- list(diag(c(1, -1, -1)), diag(c(-1, 1, -1)), diag(c(-1, -1, 1)))
  fit <- riem.mean(wrap.rotation(data), weight = c(.4, .35, .25), geometry = "extrinsic")
  # The weighted ambient matrix has negative determinant; ordinary polar
  # projection would return an improper rotation.
  expect_equal(fit$mean, data[[1]], tolerance = 1e-12)
  expect_equal(det(fit$mean), 1, tolerance = 1e-12)
  expect_error(riem.mean(wrap.rotation(data), geometry = "extrinsic"), "nonunique")
  rotate <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  angles <- c(0, .2, .5)
  fit <- riem.mean(wrap.rotation(lapply(angles, rotate)), weight = c(.2, .3, .5), eps = 1e-10)
  expect_true(fit$converged)
  expect_equal(fit$mean, rotate(sum(angles * c(.2, .3, .5))), tolerance = 1e-9)
})

test_that("Fisher simplex maps are stable at identity and on the positive domain", {
  probability <- function(a) matrix(c(cos(a)^2, sin(a)^2), 2)
  X <- probability(.3)
  Y <- probability(.9)
  identity <- geometry_operations("multinomial", X, X, 0 * X, 0 * X)
  expect_equal(identity$distance, 0)
  expect_equal(identity$log, 0 * X)
  expect_equal(identity$exp, X)
  first <- geometry_operations("multinomial", X, Y, 0 * X, 0 * X)
  maps <- geometry_operations("multinomial", X, Y, first$log, first$log)
  expect_equal(first$distance, 1.2, tolerance = 1e-12)
  expect_equal(maps$metric, first$distance^2, tolerance = 1e-12)
  expect_equal(maps$exp, Y, tolerance = 1e-11)
  expect_equal(sum(first$log), 0, tolerance = 1e-14)
  tiny <- wrap.multinomial(list(as.numeric(probability(.6)), as.numeric(probability(.6 + 1e-9))))
  expect_equal(riem.pdist(tiny)[1, 2] / 2e-9, 1, tolerance = 1e-6)
  boundary <- matrix(c(.5, .5), 2)
  expect_error(geometry_operations("multinomial", boundary, boundary,
    matrix(c(1, -1), 2), matrix(c(1, -1), 2)), "boundary")
  points <- wrap.multinomial(lapply(c(.2, .6, 1.1), function(a) as.numeric(probability(a))))
  fit <- riem.mean(points, eps = 1e-10)
  expect_true(fit$converged)
  expect_equal(fit$mean, probability(mean(c(.2, .6, 1.1))), tolerance = 1e-9)
})

test_that("fixed-rank horizontal geometry preserves the base point and equivalence class", {
  X <- matrix(c(2, 0, 0, 0, 1, 0), 3)
  Y <- matrix(c(2, 0, .3, 0, 1.5, 0), 3)
  first <- geometry_operations("spdk", X, Y, 0 * X, 0 * X)
  expect_equal(first$exp, X)
  expect_equal(first$distance, sqrt(.3^2 + .5^2), tolerance = 1e-12)
  maps <- geometry_operations("spdk", X, Y, first$log, first$log)
  expect_equal(tcrossprod(maps$exp), tcrossprod(Y), tolerance = 1e-11)
  expect_equal(maps$metric, first$distance^2, tolerance = 1e-12)
  Q <- matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2)
  changed <- geometry_operations("spdk", X %*% Q, Y %*% diag(c(-1, 1)), 0 * X, 0 * X)
  expect_equal(changed$distance, first$distance, tolerance = 1e-12)
  expect_equal(changed$log, first$log %*% Q, tolerance = 1e-11)
  wrapped <- wrap.spdk(list(tcrossprod(X), tcrossprod(Y)), k = 2)
  fit <- riem.mean(wrapped, eps = 1e-9)
  expect_true(fit$converged)
  expect_equal(tcrossprod(fit$mean), tcrossprod((X + Y) / 2), tolerance = 1e-8)

  a <- matrix(c(1, 0), 2)
  b <- matrix(c(0, 1), 2)
  expect_equal(riem.pdist(wrap.spdk(list(tcrossprod(a), tcrossprod(b)), k = 1))[1, 2], sqrt(2))
  expect_error(geometry_operations("spdk", a, b, 0 * a, 0 * a), "nonunique")
  expect_error(geometry_operations("spdk", a, a, -2 * a, -2 * a), "local horizontal domain")
})

test_that("regular landmark operations respect the orthogonal quotient", {
  raw <- rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  changed <- raw
  changed[1, ] <- changed[1, ] + c(.2, .1)
  data <- wrap.landmark(list(raw, changed))
  X <- data$data[[1]]
  Y <- data$data[[2]]
  first <- geometry_operations("landmark", X, Y, 0 * X, 0 * X)
  maps <- geometry_operations("landmark", X, Y, first$log, first$log)
  expect_equal(tcrossprod(maps$exp), tcrossprod(Y), tolerance = 1e-10)
  expect_equal(maps$metric, first$distance^2, tolerance = 1e-11)
  Q <- matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2)
  other <- geometry_operations("landmark", X %*% Q, Y %*% diag(c(-1, 1)), 0 * X, 0 * X)
  expect_equal(other$distance, first$distance, tolerance = 1e-11)
  expect_equal(other$log, first$log %*% Q, tolerance = 1e-10)
  fit <- riem.mean(data, eps = 1e-9)
  alternate <- riem.mean(wrap.landmark(list(raw %*% Q, changed %*% diag(c(-1, 1)))), eps = 1e-9)
  expect_true(fit$converged)
  expect_equal(tcrossprod(fit$mean), tcrossprod(alternate$mean), tolerance = 1e-8)

  H <- qr.Q(qr(stats::contr.helmert(5)))
  A <- H[, 1:2] / sqrt(2)
  B <- H[, 3:4] / sqrt(2)
  expect_equal(riem.pdist(wrap.landmark(list(A, B)))[1, 2], pi / 2, tolerance = 1e-12)
  expect_error(geometry_operations("landmark", A, B, 0 * A, 0 * A), "nonunique")
})

test_that("local intrinsic medians agree with independent weighted line medians", {
  weight <- c(.4, .3, .3)
  angle <- c(.1, .3, .6)
  rotate <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  line <- function(a) matrix(c(cos(a), sin(a), 0), 3)
  probability <- function(a) c(cos(a)^2, sin(a)^2)
  cases <- list(
    list(data = wrap.rotation(lapply(angle, rotate)), target = rotate(angle[2]), quotient = FALSE),
    list(data = wrap.grassmann(lapply(angle, line)), target = line(angle[2]), quotient = TRUE),
    list(data = wrap.multinomial(lapply(angle, probability)),
         target = matrix(probability(angle[2]), 2), quotient = FALSE),
    list(data = wrap.spdk(lapply(c(1, 1.5, 2), function(a) diag(c(a^2, 0, 0))), k = 1),
         target = matrix(c(1.5, 0, 0), 3), quotient = TRUE)
  )
  for (case in cases) {
    fit <- riem.median(case$data, weight = weight, eps = 1e-9)
    expect_true(fit$converged, info = case$data$name)
    observed <- if (case$quotient) tcrossprod(fit$median) else fit$median
    expected <- if (case$quotient) tcrossprod(case$target) else case$target
    expect_equal(observed, expected, tolerance = 1e-8, info = case$data$name)
    expect_lte(fit$subgradient_residual, 1e-9)
  }
})

test_that("ambient orthogonal actions preserve the declared distances", {
  X <- diag(4)[, 1:2]
  Y <- rbind(diag(cos(c(.2, .35))), diag(sin(c(.2, .35))))
  ambient <- qr.Q(qr(matrix(c(1, 3, 2, 4, 2, -1, 4, 1,
    -3, 2, 1, 2, 1, 0, -2, 3), 4)))
  for (geometry in c("intrinsic", "extrinsic")) {
    original <- riem.pdist(wrap.grassmann(list(X, Y)), geometry = geometry)
    transformed <- riem.pdist(wrap.grassmann(list(ambient %*% X, ambient %*% Y)),
      geometry = geometry)
    expect_equal(transformed, original, tolerance = 1e-12)
  }

  rx <- function(a) rbind(c(1, 0, 0), c(0, cos(a), -sin(a)), c(0, sin(a), cos(a)))
  rz <- function(a) rbind(c(cos(a), -sin(a), 0), c(sin(a), cos(a), 0), c(0, 0, 1))
  rotations <- list(rx(.1), rz(.4), rx(-.2) %*% rz(.3))
  transformed <- lapply(rotations, function(x) rz(.7) %*% x %*% rx(-.5))
  for (geometry in c("intrinsic", "extrinsic")) {
    expect_equal(riem.pdist(wrap.rotation(transformed), geometry = geometry),
      riem.pdist(wrap.rotation(rotations), geometry = geometry), tolerance = 1e-12)
  }
})

test_that("Stiefel extrinsic median reports its ambient estimand", {
  X <- diag(3)[, 1:2]
  Y <- rbind(matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2), c(0, 0))
  fit <- riem.median(wrap.stiefel(list(X, X, Y)), weight = c(.2, .35, .45),
    geometry = "extrinsic", eps = 1e-10)
  expect_true(fit$converged)
  expect_equal(fit$median, X, tolerance = 1e-12)
  expect_equal(fit$estimand, "projected_ambient_geometric_median")
  expect_equal(fit$diagnostic_scope, "ambient_embedding")
  expect_true(is.na(fit$subgradient_residual))
})
