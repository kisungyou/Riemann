special_distance <- function(x, y, geometry) {
  spd.pdist(wrap.spd(list(x, y)), geometry)[1L, 2L]
}

test_that("Stein distances preserve units and finite diagonal formulas", {
  expected <- sqrt(100 * log(1.5 / sqrt(2)))
  for (s in c(1e-200, .001, 1, 100, 1e200)) {
    expect_equal(special_distance(s * diag(100), 2 * s * diag(100), "stein"),
                 expected, tolerance = 1e-12)
  }
  x <- matrix(c(2, .3, .3, 1), 2)
  y <- matrix(c(1, -.2, -.2, 3), 2)
  oracle <- sqrt(log(det((x + y) / 2)) - (log(det(x)) + log(det(y))) / 2)
  expect_equal(special_distance(x, y, "stein"), oracle, tolerance = 1e-12)
  a <- matrix(c(1, .2, -.4, 2), 2)
  congruence <- function(z) a %*% z %*% t(a)
  expect_equal(special_distance(congruence(x), congruence(y), "stein"),
               oracle, tolerance = 1e-12)
  expect_equal(special_distance(y, x, "stein"), oracle, tolerance = 1e-12)
  expect_identical(special_distance(x, x, "stein"), 0)
  # Stable scalar identity, without subtracting logarithms close to zero.
  delta <- 1e-8
  t <- log1p(delta) / 2
  close_oracle <- sqrt(log1p(2 * sinh(t / 2)^2))
  expect_equal(special_distance(diag(1), matrix(1 + delta), "stein") / close_oracle,
               1, tolerance = 1e-7)
  extreme_oracle <- sqrt(0.5 * log(1e300) - log(2))
  expect_equal(special_distance(matrix(1e-150), matrix(1e150), "stein"),
               extreme_oracle, tolerance = 1e-12)
})

test_that("Wasserstein distances use a stable factor residual", {
  x <- matrix(c(1, .2, .2, 2), 2)
  expect_identical(special_distance(x, x, "wasserstein"), 0)
  a <- c(1, 2, 3)
  b <- a + 1e-8
  # (sqrt(b)-sqrt(a)) = (b-a)/(sqrt(b)+sqrt(a)).
  close_oracle <- sqrt(sum(((b - a) / (sqrt(b) + sqrt(a)))^2))
  for (s in c(1e-200, 1, 1e200)) {
    value <- special_distance(s * diag(a), s * diag(b), "wasserstein")
    expect_gt(value, 0)
    expect_equal(value / (sqrt(s) * close_oracle), 1, tolerance = 2e-7)
  }
  y <- matrix(c(3, -.4, -.4, 1), 2)
  # The eigenvalue sum of a positive 2 x 2 square root has this closed form.
  oracle <- sqrt(sum(diag(x)) + sum(diag(y)) -
    2 * sqrt(sum(diag(x %*% y)) + 2 * sqrt(det(x) * det(y))))
  expect_equal(special_distance(x, y, "wasserstein"), oracle, tolerance = 1e-12)
  expect_equal(special_distance(y, x, "wasserstein"), oracle, tolerance = 1e-12)
  angle <- .37
  q <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  expect_equal(special_distance(q %*% x %*% t(q), q %*% y %*% t(q), "wasserstein"),
               oracle, tolerance = 1e-12)
  for (s in c(1e-200, 1e200)) {
    expect_equal(special_distance(s * x, s * y, "wasserstein") / sqrt(s),
                 oracle, tolerance = 1e-12)
  }
})

test_that("special SPD distance output options are validated", {
  x <- wrap.spd(list(diag(2), 2 * diag(2)))
  for (geometry in c("stein", "wasserstein", "airm", "lerm")) {
    expect_s3_class(spd.pdist(x, geometry, as.dist = TRUE), "dist")
    for (bad in list(NA, 1, c(TRUE, FALSE), logical())) {
      expect_error(spd.pdist(x, geometry, as.dist = bad), "as.dist")
    }
  }
})
