test_that("weighted Euclidean means and diagonal SPD means match closed forms", {
  z <- rbind(c(-2, 4), c(1, 2), c(5, -1))
  w <- c(2, 3, 5)
  fit <- riem.mean(wrap.euclidean(z), weight = w, maxiter = 1)
  expect_s3_class(fit, "riem_summary")
  expect_true(fit$converged)
  expect_identical(fit$termination, "closed_form")
  expect_equal(as.numeric(fit$mean), colSums(z * w) / sum(w))
  expect_equal(fit$weights, w / sum(w))

  values <- rbind(c(1, 2), c(4, 8), c(9, 3))
  x <- wrap.spd(lapply(seq_len(nrow(values)), function(i) diag(values[i, ])))
  target <- diag(exp(colSums(log(values) * w) / sum(w)))
  for (geometry in c("affine_invariant", "log_euclidean")) {
    fit <- riem.mean(x, weight = w, geometry = geometry, eps = 1e-10)
    expect_true(fit$converged)
    expect_equal(unname(fit$mean), target, tolerance = 1e-8)
    expect_identical(fit$geometry$geometry_id, geometry)
    expected <- sum((w / sum(w)) * rowSums((log(values) -
      matrix(log(diag(target)), nrow(values), 2, byrow = TRUE))^2))
    expect_equal(fit$objective, expected, tolerance = 1e-10)
    expect_identical(fit$objective, fit$variation)
  }
})

test_that("affine-invariant means agree with an independently optimized reference", {
  data <- list(diag(c(1, 4)), matrix(c(3, 1, 1, 1), 2),
               matrix(c(2, -.5, -.5, 5), 2))
  w <- c(.2, .3, .5)
  fit <- riem.mean(wrap.spd(data), weight = w, eps = 1e-9, maxiter = 200, trace = TRUE)
  # The independent reference optimizes log-Cholesky coordinates using only
  # base R. It does not call a Riemann distance, logarithm, or mean routine.
  from_coordinates <- function(a) {
    L <- matrix(c(exp(a[1]), a[2], 0, exp(a[3])), 2)
    tcrossprod(L)
  }
  cost <- function(a) {
    X <- from_coordinates(a)
    sum(w * vapply(data, function(Y) {
      values <- eigen(solve(X, Y), only.values = TRUE)$values
      sum(log(values)^2)
    }, numeric(1)))
  }
  reference <- stats::optim(c(0, 0, 0), cost, method = "BFGS",
                            control = list(reltol = 1e-13, maxit = 1000))
  expect_equal(reference$convergence, 0L)
  expect_true(fit$converged)
  expect_equal(unname(fit$mean), from_coordinates(reference$par), tolerance = 2e-6)
  expect_equal(fit$objective, reference$value, tolerance = 1e-9)
  expect_lte(fit$gradient_norm, 1e-9)
  expect_true(all(diff(fit$trace$objective) <= 1e-12))
  expect_equal(tail(fit$trace$objective, 1), fit$objective)
  expect_equal(nrow(fit$trace), fit$iterations + 1L)
})

test_that("summary weights preserve duplication, order, and zero-weight exclusion", {
  a <- diag(c(1, 2))
  b <- matrix(c(4, 1, 1, 3), 2)
  base <- riem.mean(wrap.spd(list(a, b)), weight = c(.4, .6), eps = 1e-9)
  duplicate <- riem.mean(wrap.spd(list(b, a, a)), weight = c(.6, .1, .3), eps = 1e-9)
  expect_equal(base$mean, duplicate$mean, tolerance = 1e-8)
  expect_equal(base$objective, duplicate$objective, tolerance = 1e-10)

  # The zero-weight antipode must never be sent through a sphere logarithm.
  x <- wrap.sphere(rbind(c(1, 0), c(-1, 0)))
  for (statistic in c("mean", "median")) {
    fun <- get(paste0("riem.", statistic))
    fit <- fun(x, weight = c(1, 0))
    expect_true(fit$converged)
    expect_equal(as.numeric(fit[[statistic]]), c(1, 0))
    expect_equal(fit$objective, 0)
  }
  huge <- riem.mean(wrap.euclidean(rbind(c(0, 1), c(2, 3))),
                    weight = c(1e308, 1e308))
  expect_equal(as.numeric(huge$mean), c(1, 2))
})

test_that("iteration limits report an unfinished noncommuting problem", {
  data <- list(diag(c(1, 4)), matrix(c(3, 1, 1, 1), 2),
               matrix(c(2, -.5, -.5, 5), 2))
  expect_warning(fit <- riem.mean(wrap.spd(data), maxiter = 1, eps = 1e-12,
                                  trace = TRUE), "did not converge")
  expect_false(fit$converged)
  expect_identical(fit$termination, "maxiter")
  expect_equal(fit$iterations, 1L)
  expect_true(all(is.finite(fit$mean)))
  expect_gt(fit$gradient_norm, 1e-12)
  expect_equal(tail(fit$trace$objective, 1), fit$objective)
})

test_that("a bounded failed line search retains the last accepted SPD estimate", {
  rotation <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  data <- lapply(c(0, .7, 1.5), function(a) {
    Q <- rotation(a)
    Q %*% diag(c(100, .01)) %*% t(Q)
  })
  x <- wrap.spd(data)
  expect_warning(fit <- riem.mean(x, weight = c(.2, .3, .5),
    max_backtrack = 1, maxiter = 100, eps = 1e-10, trace = TRUE), "did not converge")
  expect_false(fit$converged)
  expect_identical(fit$termination, "line_search_failed")
  expect_equal(fit$iterations, 1L)
  expect_true(all(eigen(fit$mean, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(tail(fit$trace$objective, 1), fit$objective)
  expect_true(all(diff(fit$trace$objective) <= 0))
})

test_that("median coincidences retain weight and use a nonsmooth certificate", {
  x <- wrap.euclidean(rbind(c(-2, 0), c(0, 0), c(1, 0)))
  # No individual weight has a majority. The weighted initializer equals the
  # middle observation, whose coincident mass is essential to stationarity.
  fit <- riem.median(x, weight = c(.2, .4, .4), eps = 1e-12, trace = TRUE)
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$median), c(0, 0))
  expect_equal(fit$objective, .8)
  expect_lte(fit$subgradient_residual, 1e-12)

  duplicate <- wrap.euclidean(rbind(c(-2, 0), c(0, 0), c(0, 0), c(1, 0)))
  split <- riem.median(duplicate, weight = c(.2, .15, .25, .4), eps = 1e-12)
  expect_equal(split$median, fit$median)
  expect_equal(split$objective, fit$objective)

  original <- wrap.euclidean(rbind(c(-1, 0), c(0, 0), c(2, 0)))
  for (geometry in c("intrinsic", "extrinsic")) {
    answer <- riem.median(original, weight = c(.2, .7, .1), geometry = geometry,
                          maxiter = 500, eps = 1e-12)
    expect_true(answer$converged)
    expect_equal(as.numeric(answer$median), c(0, 0))
    expect_equal(answer$objective, .4)
  }
})

test_that("modified Weiszfeld leaves a nonoptimal observation and descends", {
  points <- rbind(c(0, 0), c(2, 0), c(1, sqrt(3)))
  fit <- riem.median(wrap.euclidean(points), init = matrix(c(0, 0), 2),
                     eps = 1e-9, maxiter = 300, trace = TRUE)
  expect_true(fit$converged)
  expect_gt(fit$iterations, 0)
  expect_equal(as.numeric(fit$median), c(1, 1 / sqrt(3)), tolerance = 1e-7)
  expect_lte(fit$subgradient_residual, 1e-9)
  expect_true(all(diff(fit$trace$objective) <= 1e-12))
  expect_equal(fit$objective, 2 / sqrt(3), tolerance = 1e-9)
})

test_that("median objective acceptance respects the scale of the observations", {
  triangle <- rbind(c(0, 0), c(2, 0), c(1, sqrt(3)))
  for (scale in c(1e-15, 1, 1e15)) {
    fit <- riem.median(wrap.euclidean(scale * triangle), init = matrix(c(0, 0), 2),
      eps = 1e-9, maxiter = 300, trace = TRUE)
    expect_true(fit$converged)
    expect_equal(as.numeric(fit$median) / scale, c(1, 1 / sqrt(3)), tolerance = 1e-7)
    expect_true(all(diff(fit$trace$objective) / scale <= 1e-12))
    expect_equal(fit$objective / scale, 2 / sqrt(3), tolerance = 1e-9)
  }
})

test_that("SPD median chart and intrinsic diagonal solutions are consistent", {
  x <- wrap.spd(lapply(c(-1, 0, 2), function(a) exp(a) * diag(2)))
  for (geometry in c("affine_invariant", "log_euclidean")) {
    fit <- riem.median(x, geometry = geometry, weight = c(.4, .3, .3), eps = 1e-9)
    expect_true(fit$converged)
    expect_equal(unname(fit$median), diag(2), tolerance = 1e-8)
    expect_equal(fit$objective, sqrt(2), tolerance = 1e-8)
  }
})

test_that("Grassmann projected medians use the embedding size and name their estimand", {
  basis <- function(a) cbind(c(1, 0, 0), c(0, cos(a), sin(a)))
  data <- lapply(c(0, .2, .5), basis)
  x <- wrap.grassmann(data)
  fit <- riem.median(x, geometry = "extrinsic", eps = 1e-9, maxiter = 300)
  expect_true(fit$converged)
  expect_equal(dim(fit$median), c(3L, 2L))
  expect_equal(crossprod(fit$median), diag(2), tolerance = 1e-10)
  expect_identical(fit$estimand, "projected_ambient_geometric_median")
  expect_identical(fit$diagnostic_scope, "ambient_embedding")
  expect_true(is.finite(fit$ambient_objective))
  expect_true(is.na(fit$subgradient_residual))
  rotated <- wrap.grassmann(lapply(data, function(X) X %*% diag(c(-1, 1))))
  other <- riem.median(rotated, geometry = "extrinsic", eps = 1e-9, maxiter = 300)
  expect_equal(tcrossprod(fit$median), tcrossprod(other$median), tolerance = 1e-8)
})

test_that("controls are validated and saved summary geometry is persistent", {
  x <- wrap.spd(list(diag(2), 4 * diag(2)))
  for (bad in c(0, -1, .5, Inf, NA_real_)) {
    expect_error(riem.mean(x, maxiter = bad), "positive integer")
  }
  expect_error(riem.mean(x, eps = 0), "eps")
  expect_error(riem.median(x, max_backtrack = 0), "max_backtrack")
  expect_error(riem.mean(x, made_up = 1), "Unknown summary control")
  expect_error(riem.mean(x, trace = NA), "trace")
  expect_error(riem.mean(x, init = matrix(0, 3, 3)), "init")
  expect_error(riem.mean(x, weight = c(0, 0)), "weight")
  expect_error(riem.mean(x, weight = c(1, Inf)), "weight")
  fit <- riem.mean(x, geometry = "log_euclidean")
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(fit, file)
  expect_identical(readRDS(file), fit)
  expect_identical(fit$geometry$geometry_id, "log_euclidean")
  expect_identical(fit$schema_version, 1L)
})
