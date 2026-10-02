test_that("legacy wrappers reject invalid membership before normalization", {
  expect_error(wrap.multinomial(list(c(-1, -2))), "positive")
  expect_error(wrap.multinomial(list(c(0, 2))), "positive")
  expect_error(wrap.multinomial(list(c(NA, 2))), "finite")
  expect_equal(as.numeric(wrap.multinomial(list(c(1e300, 2e300)))$data[[1]]), c(1/3, 2/3))
  expect_error(wrap.correlation(list(matrix(c(1, 2, 2, 1), 2))), "positive definite")
  expect_error(wrap.correlation(list(diag(c(2, 1)))), "unit diagonal")
  expect_error(wrap.rotation(list(diag(c(1, -1)))), "determinant")
  expect_error(wrap.stiefel(list(matrix(1, 3, 2))), "full-column-rank")
  expect_error(wrap.grassmann(list(matrix(1, 2, 3))), "full-column-rank")
  expect_equal(wrap.grassmann(array(c(1, 0, 0), c(3, 1, 1)))$size, c(3L, 1L))
  expect_equal(wrap.rotation(array(1, c(1, 1, 1)))$size, c(1L, 1L))
})

test_that("fixed-rank wrapping requires positive semidefiniteness and exact numerical rank", {
  x <- diag(c(3, 1, 0))
  z <- wrap.spdk(list(x), 2)$data[[1]]
  expect_equal(tcrossprod(z), x, tolerance = 1e-12)
  z1 <- wrap.spdk(list(diag(c(3, 0, 0))), 1)$data[[1]]
  expect_equal(dim(z1), c(3L, 1L))
  expect_error(wrap.spdk(list(diag(3)), 2), "rank exactly")
  expect_error(wrap.spdk(list(diag(c(3, 1, -1))), 2), "positive semidefinite")
  expect_error(wrap.spdk(list(x), 1.5), "integer")
  expect_error(wrap.spdk(list(x), 3), "p-1")
  zsmall <- wrap.spdk(list(x * 1e-200), 2)$data[[1]]
  expect_equal(tcrossprod(zsmall / 1e-100), x, tolerance = 1e-12)
})

test_that("unavailable legacy geometries are rejected explicitly", {
  st <- wrap.stiefel(list(diag(3)[, 1:2]))
  expect_error(riem.pdist(st), "does not support distance")
  expect_error(riem.mean(st), "does not support mean")
  corr <- wrap.correlation(list(diag(2)))
  expect_error(riem.pdist(corr), "does not support distance")
  expect_error(riem.mean(corr), "does not support mean")
  psd <- wrap.spdk(list(diag(c(2, 1, 0))), 2)
  expect_error(riem.pdist(psd, geometry = "extrinsic"), "does not support distance")
})

test_that("saved models reject contradictory geometry instead of changing interpretation", {
  x <- wrap.spd(list(diag(c(1, 2)), diag(c(2, 1)), diag(c(3, 2))))
  pca <- riem.pga(x, ndim = 1, geometry = "log_euclidean")
  cluster <- riem.kmeans(x, k = 2, geometry = "log_euclidean", init = c(1, 2), nstart = 1)
  reg <- riem.m2skreg(x, 1:3, geometry = "log_euclidean")
  for (model in list(pca, cluster, reg)) {
    file <- tempfile(fileext = ".rds")
    saveRDS(model, file)
    expect_equal(predict(readRDS(file), x), predict(model, x))
    bad <- model
    bad$geometry$backend <- "intrinsic"
    expect_error(predict(bad, x), "metadata")
    bad <- model
    bad$geometry$schema_version <- 99L
    expect_error(predict(bad, x), "schema")
    unlink(file)
  }
  pca$coordinates$geometry <- riem.geometry(x, "intrinsic")
  expect_error(predict(pca, x), "incompatible")
})

test_that("trained summaries preserve named coordinates through model reuse", {
  expect_error(wrap.euclidean(list(c(1, 100), c(b = 100, a = 1))), "consistent feature")
  expect_error(wrap.sphere(list(c(1, 0), c(b = 0, a = 1))), "consistent feature")
  named <- diag(2)
  dimnames(named) <- list(c("a", "b"), c("a", "b"))
  expect_error(wrap.spd(list(diag(2), named)), "consistent feature")
  x <- wrap.euclidean(matrix(c(1, 3, 6, 2, 0, 4), ncol = 2,
                              dimnames = list(NULL, c("length", "width"))))
  expect_equal(rownames(riem.mean(x)$mean), c("length", "width"))
  pca <- riem.pga(x, ndim = 1)
  expect_equal(predict(pca, x), pca$embed)
  clustering <- riem.kmeans(x, k = 2, init = c(1, 3), nstart = 1)
  expect_equal(predict(clustering, x), clustering$cluster)
  swapped <- wrap.euclidean(matrix(c(2, 0, 4, 1, 3, 6), ncol = 2,
                                    dimnames = list(NULL, c("width", "length"))))
  expect_error(predict(pca, swapped), "incompatible")
  expect_error(predict(clustering, swapped), "incompatible")
})
