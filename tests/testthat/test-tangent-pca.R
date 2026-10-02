test_that("tangent PCA agrees with ordinary PCA and preserves requested dimensions", {
  z <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 2, 0),
             c(0, 0, 3), c(1, 1, 1), c(-1, 2, -2))
  x <- wrap.euclidean(z)
  reference <- stats::prcomp(z)
  for (k in 1:3) {
    fit <- riem.pga(x, ndim = k)
    expect_equal(dim(fit$embed), c(nrow(z), k))
    expect_equal(fit$variance, reference$sdev[seq_len(k)]^2, tolerance = 1e-9)
    expect_equal(tcrossprod(fit$embed), tcrossprod(reference$x[, seq_len(k), drop = FALSE]),
                 tolerance = 1e-9)
    expect_equal(unname(predict(fit, x)), unname(fit$embed), tolerance = 1e-12)
  }
  reconstructed <- riem.reconstruct(fit)
  expect_equal(do.call(rbind, lapply(reconstructed$data, as.numeric)), z, tolerance = 1e-10)
  expect_true(fit$diagnostics$converged)
  expect_null(fit$input_template$data)
})

test_that("weighted tangent PCA matches an independently assembled covariance", {
  z <- rbind(c(1, 2), c(3, -1), c(-2, 1), c(1000, 1000))
  w <- c(0.2, 0.5, 0.3, 0)
  fit <- riem.pga(wrap.euclidean(z), ndim = 2, weight = w)
  mu <- colSums(z * w)
  centered <- sweep(z, 2, mu)
  covariance <- Reduce(`+`, lapply(seq_len(nrow(z)), function(i) {
    w[i] * tcrossprod(centered[i, ])
  })) / (1 - sum(w^2))
  expected <- eigen(covariance, symmetric = TRUE)$values
  expect_equal(as.numeric(fit$center), mu, tolerance = 1e-10)
  expect_equal(fit$variance, expected, tolerance = 1e-10)
  unweighted.outlier <- riem.pga(wrap.euclidean(z[1:3, ]), ndim = 2, weight = w[1:3])
  expect_equal(fit$variance, unweighted.outlier$variance, tolerance = 1e-10)
  expect_equal(fit$loadings, unweighted.outlier$loadings, tolerance = 1e-10)
})

test_that("symmetric coordinates preserve Frobenius and affine-invariant metrics", {
  u <- matrix(c(2, 3, 3, -1), 2)
  v <- matrix(c(1, -2, -2, 4), 2)
  expect_equal(riem_tangent_svec(u), c(2, sqrt(2) * 3, -1))
  expect_equal(sum(riem_tangent_svec(u) * riem_tangent_svec(v)), sum(u * v))
  c <- matrix(c(2, 0.4, 0.4, 1), 2)
  invroot <- riem_tangent_symmetric_function(c, function(x) 1 / sqrt(x), positive = TRUE)
  whitened.u <- invroot %*% u %*% invroot
  whitened.v <- invroot %*% v %*% invroot
  expected <- sum(diag(solve(c, u) %*% solve(c, v)))
  expect_equal(sum(riem_tangent_svec(whitened.u) * riem_tangent_svec(whitened.v)),
               expected, tolerance = 1e-12)
})

test_that("affine-invariant diagonal PGA has the independent log-coordinate answer", {
  logs <- rbind(c(0, 0), c(0.5, 1.2), c(-0.4, 0.9), c(0.8, -0.1), c(-0.6, -0.5))
  x <- wrap.spd(lapply(seq_len(nrow(logs)), function(i) diag(exp(logs[i, ]))))
  fit <- riem.pga(x, ndim = 2, geometry = "affine_invariant", eps = 1e-10)
  expected <- stats::prcomp(logs)
  expect_equal(fit$center, diag(exp(colMeans(logs))), tolerance = 1e-8)
  expect_equal(tcrossprod(fit$embed), tcrossprod(expected$x), tolerance = 1e-8)
  expect_equal(riem.reconstruct(fit)$data, x$data, tolerance = 1e-8)
  expect_identical(fit$geometry$geometry_id, "affine_invariant")
})

test_that("log-Euclidean PCA handles noncommuting matrices with chart-correct covariance", {
  charts <- list(matrix(c(0.2, 0.1, 0.1, -0.3), 2),
                 matrix(c(0.9, -0.3, -0.3, 0.4), 2),
                 matrix(c(-0.4, 0.2, 0.2, 0.8), 2),
                 matrix(c(0.1, 0.6, 0.6, 0.4), 2),
                 matrix(c(0.3, -0.4, -0.4, -0.7), 2))
  # Independent spectral construction, not the production matrix-function helper.
  exponential <- function(a) {
    e <- eigen(a, symmetric = TRUE)
    e$vectors %*% diag(exp(e$values)) %*% t(e$vectors)
  }
  x <- wrap.spd(lapply(charts, exponential))
  fit <- riem.pga(x, ndim = 3, geometry = "log_euclidean")
  expected <- stats::prcomp(do.call(rbind, lapply(charts, function(a) {
    c(a[1, 1], sqrt(2) * a[1, 2], a[2, 2])
  })))
  expect_equal(fit$center, exponential(Reduce(`+`, charts) / length(charts)),
               tolerance = 1e-10)
  expect_equal(tcrossprod(fit$embed), tcrossprod(expected$x), tolerance = 1e-10)
  expect_equal(fit$variance, expected$sdev^2, tolerance = 1e-10)
  expect_equal(riem.reconstruct(fit)$data, x$data, tolerance = 1e-10)
  expect_error(predict(fit, wrap.euclidean(matrix(1:6, 3))), "incompatible")
})

test_that("sphere tangent representation is metric-correct and prediction is batch invariant", {
  z <- rbind(c(1, 0, 0), c(1, 0.2, 0), c(1, -0.2, 0),
             c(1, 0, 0.3), c(1, 0, -0.3))
  x <- wrap.sphere(z)
  fit <- riem.pga(x, ndim = 2, eps = 1e-10)
  expect_equal(crossprod(fit$coordinates$basis), diag(2), tolerance = 1e-12)
  expect_equal(as.numeric(crossprod(fit$center, fit$coordinates$basis)), c(0, 0),
               tolerance = 1e-12)
  expect_equal(riem.reconstruct(fit)$data, x$data, tolerance = 1e-9)
  batch <- predict(fit, x)
  singleton <- do.call(rbind, lapply(seq_len(nrow(z)), function(i) {
    predict(fit, wrap.sphere(z[i, , drop = FALSE]))
  }))
  expect_equal(singleton, batch, tolerance = 1e-12)
  order <- c(5, 1, 3, 2, 4)
  expect_equal(predict(fit, wrap.sphere(z[order, ])), batch[order, ], tolerance = 1e-12)
  expect_error(predict(fit, wrap.sphere(matrix(c(-1, 0, 0), 1))), "antipode")
  too.far <- matrix(c(4, 0), 1)
  expect_error(riem.reconstruct(fit, too.far), "injectivity")
})

test_that("rank truncation, constant data, persistence, and unsupported cases are explicit", {
  x <- wrap.euclidean(cbind(1:4, 2 * (1:4)))
  expect_warning(fit <- riem.pga(x, ndim = 2), "numerical rank")
  expect_equal(fit$rank, 1)
  expect_equal(dim(fit$embed), c(4, 1))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(predict(restored, x), predict(fit, x))
  expect_identical(restored$geometry, fit$geometry)
  expect_warning(constant <- riem.pga(wrap.euclidean(matrix(1, 3, 2))), "numerical rank")
  expect_equal(dim(constant$embed), c(3, 0))
  expect_equal(riem.reconstruct(constant)$data, rep(list(matrix(1, 2, 1)), 3))
  expect_error(riem.pga(x, ndim = 1.5), "positive integer")
  expect_error(riem.pga(x, rank.tol = 0), "rank.tol")
  expect_error(riem.pga(wrap.sphere(cbind(rep(1, 3), 1:3)), geometry = "chordal"),
               "does not support tangent")
  fit$schema_version <- 999L
  expect_error(predict(fit, x), "supported fitted")
})

test_that("regular landmark prediction retains the trained orthogonal-quotient frame", {
  base <- rbind(c(-1, -0.6), c(0.8, -0.5), c(0.9, 0.7), c(-0.7, 0.8))
  set.seed(9021)
  raw <- lapply(seq_len(8), function(i) base + matrix(stats::rnorm(8, sd = 0.05), 4, 2))
  x <- wrap.landmark(raw)
  fit <- riem.pga(x, ndim = 4, eps = 1e-9)
  original <- predict(fit, x)
  angle <- 0.8
  rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  transformed <- lapply(raw, function(a) sweep(3 * a %*% rotation, 2, c(4, -2), "+"))
  reflected <- lapply(raw, function(a) a %*% diag(c(-1, 1)))
  expect_equal(predict(fit, wrap.landmark(transformed)), original, tolerance = 1e-7)
  expect_equal(predict(fit, wrap.landmark(reflected)), original, tolerance = 1e-7)
  singleton <- do.call(rbind, lapply(raw, function(a) predict(fit, wrap.landmark(list(a)))))
  expect_equal(singleton, original, tolerance = 1e-7)
  reconstructed <- riem.reconstruct(fit)
  for (i in seq_along(x$data)) {
    expect_equal(colSums(reconstructed$data[[i]]), c(0, 0), tolerance = 1e-10)
    expect_equal(sum(reconstructed$data[[i]]^2), 1, tolerance = 1e-10)
    aligned <- riem_tangent_landmark_align(fit$center, x$data[[i]])
    expect_equal(reconstructed$data[[i]], aligned, tolerance = 1e-7)
  }
  expect_identical(fit$coordinates$alignment.group, "O(p)")
})
test_that("landmark prediction rejects cut loci independently of representative", {
  basis <- qr.Q(qr(stats::contr.helmert(5)))
  center <- basis[, 1:2] / sqrt(2)
  orthogonal <- basis[, 3:4] / sqrt(2)
  training <- wrap.landmark(list(cos(.1) * center + sin(.1) * orthogonal,
    cos(.1) * center - sin(.1) * orthogonal))
  fit <- riem.pga(training, ndim = 1, eps = 1e-12)
  rotation <- matrix(c(cos(.2), sin(.2), -sin(.2), cos(.2)), 2)
  equivalent <- orthogonal %*% rotation
  expect_equal(riem.pdist(wrap.landmark(list(orthogonal, equivalent)))[1, 2],
    0, tolerance = 1e-12)
  expect_error(predict(fit, wrap.landmark(list(orthogonal))), "singular|nonunique")
  expect_error(predict(fit, wrap.landmark(list(equivalent))), "singular|nonunique")
})

test_that("identical SPD observations have zero tangent rank at every repetition count", {
  observation <- matrix(c(2, .3, .3, 1), 2)
  for (n in c(3L, 5L, 7L)) {
    data <- wrap.spd(rep(list(observation), n))
    for (geometry in c("affine_invariant", "log_euclidean")) {
      expect_warning(fit <- riem.pga(data, ndim = 1, geometry = geometry),
        "numerical rank is 0")
      expect_equal(fit$rank, 0)
      expect_equal(fit$total.variance, 0)
      expect_equal(dim(fit$embed), c(n, 0L))
      expect_equal(predict(fit, data), fit$embed)
      saved <- tempfile(fileext = ".rds")
      saveRDS(fit, saved)
      restored <- readRDS(saved)
      unlink(saved)
      expect_equal(dim(predict(restored, data)), c(n, 0L))
      expect_equal(predict(restored, data), fit$embed)
      expect_equal(riem.reconstruct(restored)$data, data$data, tolerance = 1e-12)
    }
  }
})
