test_that("SPD FANOVA distinguishes equal centers from equal centers and dispersions", {
  # Along this commuting SPD geodesic, the log of the first diagonal entry is
  # an isometric scalar coordinate for both metrics. Squared distances and
  # both terms of the published statistic are therefore known independently.
  wide <- c(-1, -.5, -.2, -.1, .1, .2, .5, 1)
  narrow <- wide / 8
  make_group <- function(z) wrap.spd(lapply(z, function(t) diag(c(exp(t), 1))))
  x <- make_group(wide)
  y <- make_group(narrow)
  # Equal sample sizes and identical log-coordinate means reduce equation (14)
  # to a two-sample squared difference in dispersions with a known denominator.
  v <- mean(wide^2)
  sigma2 <- mean((wide^2 - v)^2)
  expected <- length(wide) * (v - v / 64)^2 / (sigma2 + sigma2 / 64^2)
  for (geometry in c("affine_invariant", "log_euclidean")) {
    null <- riem.fanova(x, x, geometry = geometry)
    expect_equal(null$p.value, 1, tolerance = 1e-10)
    fit <- riem.fanova(x, y, geometry = geometry)
    expect_equal(fit$terms, c(variance = expected, mean = 0), tolerance = 1e-10)
    expect_equal(fit$p.value, pchisq(expected, 1L, lower.tail = FALSE))
  }
})
