# Three separated groups with distinct observations keep the backend interfaces
# under test independent of tie-breaking or arbitrary cluster identifiers.
backend_positions <- c(-2.10, -2.03, -1.97, -1.91, -.12, -.04,
                       .03, .11, 1.87, 1.95, 2.03, 2.12)

test_that("PAM medoids and competitive quantization label the nearest returned center", {
  z <- backend_positions
  data <- wrap.euclidean(as.list(z))
  pam <- riem.kmedoids(data, k = 3)
  expect_length(unique(pam$medoids), 3)
  distances <- abs(outer(z, z[pam$medoids], "-"))
  expect_identical(pam$cluster, max.col(-distances, ties.method = "first"))

  set.seed(71)
  quantized <- riem.clrq(data, k = 3, init = "plus")
  expect_equal(dim(quantized$centers), c(1, 1, 3))
  distances <- abs(outer(z, as.numeric(quantized$centers), "-"))
  expect_equal(quantized$cluster, max.col(-distances, ties.method = "first"))
})

test_that("mean shift respects the update budget and silhouette search bound", {
  z <- backend_positions
  data <- wrap.euclidean(as.list(z))
  h <- .2
  distances <- abs(outer(z, z, "-"))
  weights <- exp(-(distances / h)^2)
  updated <- as.vector(weights %*% z) / rowSums(weights)
  fit <- riem.nmshift(data, h = h, maxiter = 1, maxk = 2)
  expect_equal(fit$distance, abs(outer(updated, updated, "-")), tolerance = 1e-12)
  expect_equal(length(unique(fit$cluster)), 2)
  narrow <- riem.nmshift(data, h = 1e-200, maxiter = 1, maxk = 2)
  expect_equal(narrow$distance, distances, tolerance = 1e-12)
  expect_error(riem.nmshift(data, h = 0), "finite positive")
  expect_error(riem.nmshift(data, maxk = 1), "maxk")
})

test_that("t-SNE and PHATE distance backends return one finite row per observation", {
  data <- wrap.euclidean(as.list(backend_positions))
  set.seed(72)
  tsne <- riem.tsne(data, ndim = 2, perplexity = 2, max_iter = 100, num_threads = 1)
  phate <- riem.phate(data, ndim = 2, nbdk = 3)
  for (fit in list(tsne, phate)) {
    expect_equal(dim(fit$embed), c(length(backend_positions), 2))
    expect_true(all(is.finite(fit$embed)))
    expect_gt(max(stats::dist(fit$embed)), 0)
  }
  expect_true(length(tsne$stress) == 1L && is.finite(tsne$stress) && tsne$stress >= 0)
})
