test_that("nearest neighbors explicitly exclude self under zero-distance ties", {
  x <- wrap.euclidean(matrix(c(0, 0, 1), ncol = 1))
  fit <- riem.knn(x, 1L)
  expect_identical(as.vector(fit$nn.idx), c(2L, 1L, 1L))
  expect_equal(as.vector(fit$nn.dists), c(0, 0, 1))
  all <- riem.knn(x, 2L)
  expect_equal(all$nn.idx, rbind(c(2L, 3L), c(1L, 3L), c(1L, 2L)))
  expect_false(any(all$nn.idx == row(all$nn.idx)))
  x <- wrap.sphere(rbind(c(1, 0), c(1, 0), c(0, 1)))
  expect_equal(as.vector(riem.knn(x, 1L, "round")$nn.dists), c(0, 0, pi / 2))
  expect_identical(riem.knn(x, 1L, riem.geometry(x)), riem.knn(x, 1L))
})

test_that("nearest neighbor count and geometry have explicit contracts", {
  x <- wrap.euclidean(matrix(1:3, ncol = 1))
  for (bad in list(0, -1, 1.5, 3, Inf, NA, c(1, 2), "1")) {
    expect_error(riem.knn(x, bad), "k.*integer")
  }
  expect_error(riem.knn(wrap.euclidean(matrix(1)), 1), "at least two")
  expect_error(riem.knn(x, 1, "round"), "Unsupported geometry")
})
