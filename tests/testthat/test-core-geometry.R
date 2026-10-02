test_that("SPD validation separates definiteness from rank and physical units", {
  for (bad in list(diag(c(-1, 1)), matrix(0, 2, 2), matrix(NA_real_, 2, 2),
                   matrix(c(1, 0, 1, 1), 2), matrix(1i, 2, 2))) {
    expect_error(wrap.spd(list(bad)), "positive-definite")
    expect_error(wrap.spd(array(bad, c(2, 2, 1))), "positive-definite")
  }
  for (scale in c(1, 1e-27, 1e100)) {
    x <- diag(c(1, 1e-8)) * scale
    expect_equal(wrap.spd(list(x))$data[[1]], x)
  }
  expect_equal(wrap.spd(array(2, c(1, 1, 1)))$size, c(1L, 1L))
  expect_error(wrap.spd(list()), "nonempty")
  expect_error(wrap.spd(array(numeric(), c(2, 2, 0))), "nonempty")
})

test_that("SPD native primitives agree with AIRM and congruence invariance", {
  x <- diag(c(1, 2))
  u <- matrix(c(0, 1, 1, 0), 2)
  y <- matrix(c(2, .4, .4, 3), 2)
  op <- Riemann:::geometry_operations("spd", x, y, u, u)
  expect_equal(op$metric, 1, tolerance = 1e-12)
  z <- Riemann:::geometry_operations("spd", x, y, op$log, op$log)
  expect_equal(z$metric, z$distance^2, tolerance = 1e-10)
  expect_equal(z$exp, y, tolerance = 1e-10)
  a <- matrix(c(1, .3, -.2, 2), 2)
  cong <- function(m) a %*% m %*% t(a)
  transformed <- Riemann:::geometry_operations("spd", cong(x), cong(y), cong(u), cong(u))
  expect_equal(transformed$metric, op$metric, tolerance = 1e-10)
  expect_equal(transformed$distance, op$distance, tolerance = 1e-10)
  scaled <- Riemann:::geometry_operations("spd", x*1e-27, y*1e-27, u*1e-27, u*1e-27)
  expect_equal(scaled$metric, op$metric, tolerance = 1e-10)
  expect_equal(scaled$distance, op$distance, tolerance = 1e-10)
})

test_that("sphere primitives preserve small angles and reject nonunique logs", {
  x <- matrix(c(1, 0), 2)
  for (theta in c(1e-9, 1e-5, pi/2, pi-1e-8)) {
    y <- matrix(c(cos(theta), sin(theta)), 2)
    u <- matrix(c(0, theta), 2)
    op <- Riemann:::geometry_operations("sphere", x, y, u, u)
    expect_equal(op$distance/theta, 1, tolerance = 1e-8)
    expect_equal(op$log, u, tolerance = 1e-8)
    expect_equal(op$exp, y, tolerance = 1e-10)
  }
  expect_error(Riemann:::geometry_operations("sphere", x, -x, x*0, x*0), "antipodal")
  d <- riem.pdist(wrap.sphere(rbind(c(1, 0), c(-1, 0))))
  expect_equal(d[1, 2], pi)
  for (bad in list(c(0, 0), c(NA, 1), c(Inf, 0), c(1i, 0))) {
    expect_error(wrap.sphere(list(bad)))
  }
  expect_error(wrap.sphere(list()), "nonempty")
  huge <- wrap.sphere(rbind(c(1e300, 1e300)))
  expect_equal(sum(huge$data[[1]]^2), 1, tolerance = 1e-14)
})

test_that("zero weights and overflow-safe normalization have explicit semantics", {
  expect_equal(Riemann:::check_weight(c(0, 2, 2), 3, "test"), c(0, .5, .5))
  expect_equal(Riemann:::check_weight(c(1e308, 1e308), 2, "test"), c(.5, .5))
  for (w in list(c(0, 0), c(-1, 2), c(NA, 2), c(Inf, 1), c("1", "2"))) {
    expect_error(Riemann:::check_weight(w, 2, "test"), "finite nonnegative")
  }
})

test_that("geometry aliases persist and incompatible specifications fail", {
  x <- wrap.spd(list(diag(2), 2*diag(2)))
  g <- riem.geometry(x, "extrinsic")
  expect_identical(g$geometry_id, "log_euclidean")
  expect_identical(g, riem.geometry(x, "log_euclidean"))
  expect_identical(g, riem.geometry(x, g))
  expect_error(riem.geometry(wrap.spd(list(diag(3))), g), "dimensions")
  expect_error(riem.geometry(x, "chordal"), "Unsupported")
  expect_error(Riemann:::riem_resolve_geometry(wrap.sphere(rbind(c(1, 0))), "chordal", "tangent"),
               "does not support")
})
test_that("cached log-Euclidean distances preserve direct matrix-log differences", {
  angle <- 0.31
  rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  matrices <- list(diag(c(1, 3)), rotation %*% diag(c(2, 5)) %*% t(rotation),
                   diag(c(1 + 1e-10, 3)))
  for (scale in c(1, 1e-24, 1e100)) {
    x <- wrap.spd(lapply(matrices, function(m) m * scale))
    logs <- lapply(x$data, function(m) {
      e <- eigen(m, symmetric = TRUE)
      tcrossprod(sweep(e$vectors, 2L, log(e$values), "*"), e$vectors)
    })
    reference <- outer(seq_along(logs), seq_along(logs), Vectorize(function(i, j) {
      sqrt(sum((logs[[i]] - logs[[j]])^2))
    }))
    expect_equal(riem.pdist(x, "log_euclidean"), reference, tolerance = 1e-11)
    single <- x
    single$data <- x$data[1L]
    expect_equal(as.numeric(riem.pdist2(single, x, "log_euclidean")),
                 reference[1L, ], tolerance = 1e-11)
    expect_gt(riem.pdist(x, "log_euclidean")[1L, 3L], 0)
  }
})
