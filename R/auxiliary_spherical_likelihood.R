# Stable radial likelihood calculations shared by the spherical families.
# eta is the natural rate: lambda/2 for normal, 1/sigma for Laplace.
# The angular density is exp(-eta*r^power) sin(r)^(D-2), 0 <= r <= pi.

#' @keywords internal
#' @noRd
sphere_radial_stats <- function(eta, D, power, moments = TRUE) {
  stopifnot(length(eta) == 1L, is.finite(eta), eta >= 0,
            D >= 2, power %in% c(1, 2))
  log.scale <- if (eta > 1) log(eta)/power else 0
  scale <- exp(log.scale)
  alpha <- if (eta > 1) 1 else eta
  exponent <- D - 2
  # Scaling the angle prevents quadrature from missing a concentrated peak.
  log.kernel <- function(u) {
    z <- u/scale
    ans <- -alpha * u^power
    if (exponent > 0) {
      small <- abs(z) < 1e-4
      log.sinc <- numeric(length(z))
      log.sinc[small] <- log1p(-z[small]^2/6 + z[small]^4/120)
      log.sinc[!small] <- log(pmax(0, sin(z[!small]))/z[!small])
      ans <- ans + exponent * (log(u) + log.sinc)
    }
    ans[z > pi] <- -Inf
    ans
  }
  endpoint <- pi * scale
  if (exponent == 0) {
    mode <- 0
  } else if (alpha == 0) {
    mode <- pi/2
  } else {
    # sin(r)/r decreases on (0, pi); the gamma-kernel mode is an upper bound.
    bound <- min(endpoint/2, (exponent/(power*alpha))^(1/power))
    mode <- stats::optimize(function(u) -log.kernel(u), c(0, bound),
                            tol = 1e-10 * max(1, bound))$minimum
  }
  shift <- log.kernel(mode)
  upper <- endpoint
  if (scale > 1) {
    upper <- min(endpoint, max(1, mode + 1))
    # The log kernel is concave. Stop only after its declining tail is
    # negligible relative to its mode; the 50-log-unit margin also covers
    # the first three moments used below.
    while (upper < endpoint && log.kernel(upper) > shift - 50 - 3*power*log(max(1, upper))) {
      upper <- min(endpoint, upper * 2)
    }
  }
  weight <- function(u) exp(log.kernel(u) - shift)
  integrate.scaled <- function(f) {
    left <- if (mode > 0) stats::integrate(f, 0, mode,
      rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)$value else 0
    right <- stats::integrate(f, mode, upper,
      rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)$value
    left + right
  }
  integral <- integrate.scaled(weight)
  log.area <- log(2) + ((D-1)/2)*log(pi) - lgamma((D-1)/2)
  logZ <- log.area + shift + log(integral) - (D-1)*log.scale
  if (eta == 0) logZ <- log(2) + (D/2)*log(pi) - lgamma(D/2)
  if (!moments) return(list(logZ = logZ))
  mean.u <- integrate.scaled(function(u) u^power * weight(u))/integral
  # Central moments avoid subtracting nearly equal raw moments.
  var.u <- integrate.scaled(function(u) (u^power - mean.u)^2 * weight(u))/integral
  # The third moment can cancel, so integrate its positive and negative parts.
  third.pos <- integrate.scaled(function(u) pmax(0, u^power - mean.u)^3 * weight(u))/integral
  third.neg <- integrate.scaled(function(u) pmax(0, mean.u - u^power)^3 * weight(u))/integral
  list(logZ = logZ, mean = mean.u * exp(-power*log.scale),
       scaled.mean = alpha * mean.u, scaled.var = alpha^2 * var.u,
       scaled.third = alpha^3 * (third.pos - third.neg))
}

#' @keywords internal
#' @noRd
sphere_radial_mle <- function(target, D, power, method, maxiter, tol,
                              exact = TRUE) {
  if (target == 0) return(Inf) # point-mass limit; no finite MLE
  uniform <- sphere_radial_stats(0, D, power)$mean
  if (target >= uniform - 64*.Machine$double.eps*uniform) return(0)
  tol <- max(tol, 1e-10)
  evaluate <- function(t) sphere_radial_stats(exp(t), D, power)
  score <- function(t) exp(t)*target - evaluate(t)$scaled.mean
  # The expected statistic decreases strictly with eta (its derivative is
  # minus its variance). Hence score signs provide a valid global bracket.
  initial <- log((D-1)/power) - log(target)
  # Normal concentration is 2*eta, whereas the Laplace rate is eta itself.
  max.rate <- .Machine$double.xmax / if (power == 2) 2 else 1
  limit <- log(max.rate)
  lo <- hi <- min(initial, limit)
  while (score(lo) >= 0) lo <- lo - log(2)
  while (score(hi) < 0) {
    if (hi >= limit) {
      boundary <- evaluate(hi)
      if (abs(exp(hi)*target-boundary$scaled.mean) <= tol*boundary$scaled.mean)
        return(min(exp(hi), max.rate))
      warning("The finite radial MLE exceeds the representable parameter range.", call. = FALSE)
      return(Inf)
    }
    hi <- min(limit, hi + log(2))
  }
  loss <- function(t) exp(t)*target + sphere_radial_stats(exp(t), D, power,
                                                       moments = FALSE)$logZ
  if (method == "optimize") {
    t <- stats::optimize(loss, c(lo, hi), tol = tol)$minimum
  } else if (method == "de") {
    fit <- DEoptim::DEoptim(loss, lower = lo, upper = hi,
      control = DEoptim::DEoptim.control(trace = FALSE, itermax = maxiter,
        reltol = tol, steptol = maxiter, NP = 40L))
    t <- as.double(fit$optim$bestmem)
  } else {
    t <- (lo + hi)/2
  }
  # Newton/Halley solve the monotone score. Optimize and DE retain their
  # distinct likelihood searches, followed by the same score check and,
  # only when needed, safeguarded Newton polishing.
  converged <- FALSE
  for (iteration in seq_len(maxiter)) {
    current <- evaluate(t)
    residual <- exp(t)*target - current$scaled.mean
    if (abs(residual) <= tol * max(current$scaled.mean, exp(t)*target)) {
      converged <- TRUE
      break
    }
    if (residual < 0) lo <- t else hi <- t
    if (hi - lo <= tol) {
      t <- (lo + hi)/2
      converged <- TRUE
      break
    }
    variance <- current$scaled.var
    if (!exact && method == "newton") {
      # Finite-difference the expected statistic in log-rate coordinates.
      # This avoids the old fixed angular grid that missed narrow peaks.
      h <- 1e-4
      m.minus <- evaluate(t-h)$scaled.mean
      m.plus <- evaluate(t+h)$scaled.mean
      variance <- (exp(h)*m.minus - exp(-h)*m.plus)/(2*h)
    }
    relative.step <- -residual/variance
    if (method == "halley") {
      relative.step <- -2*residual*variance /
        (2*variance^2 + residual*current$scaled.third)
    }
    candidate <- if (is.finite(relative.step) && relative.step > -1)
      t + log1p(relative.step) else NA_real_
    if (!is.finite(candidate) || candidate <= lo || candidate >= hi) {
      candidate <- (lo + hi)/2
    }
    t <- candidate
  }
  if (!converged) warning("Radial likelihood iteration limit reached before the score tolerance.",
                           call. = FALSE)
  # exp(log(max.rate)) may round just above max.rate, which would make
  # 2*eta overflow in the normal family even though the limit is representable.
  min(exp(t), max.rate)
}

#' @keywords internal
#' @noRd
sphere_distribution_distances <- function(mu, x) {
  if (length(mu) != ncol(x) || length(mu) < 2L) {
    stop("The location and observations must have the same dimension, at least two.", call. = FALSE)
  }
  mu <- mu/sqrt(sum(mu^2))
  row.norm <- function(z) {
    size <- apply(abs(z), 1L, max)
    size.safe <- ifelse(size == 0, 1, size)
    size * sqrt(rowSums((z/size.safe)^2))
  }
  mu.rows <- matrix(mu, nrow(x), ncol(x), byrow = TRUE)
  # Half-angle chord formula: accurate at coincidence and the antipode,
  # without an acos cancellation or a tangent projection residual at x=mu.
  2*atan2(row.norm(x-mu.rows), row.norm(x+mu.rows))
}

#' @keywords internal
#' @noRd
sphere_distribution_controls <- function(pars) {
  maxiter <- if (is.null(pars$maxiter)) 50 else pars$maxiter
  eps <- if (is.null(pars$eps)) 1e-6 else pars$eps
  if (!is.numeric(maxiter) || length(maxiter) != 1L || !is.finite(maxiter) || maxiter < 1)
    stop("'maxiter' must be a positive finite number.", call. = FALSE)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    stop("'eps' must be a positive finite number.", call. = FALSE)
  list(maxiter = max(10L, round(maxiter)), eps = min(eps, 1e-6))
}

#' @keywords internal
#' @noRd
sphere_coincident_rows <- function(x) {
  all(x == matrix(x[1L, ], nrow(x), ncol(x), byrow = TRUE))
}


#' @keywords internal
#' @noRd
sphere_distribution_location <- function(spobj, x, maxiter, eps, median = FALSE) {
  if (sphere_coincident_rows(x)) return(x[1L, ])
  weights <- rep(1/nrow(x), nrow(x))
  fit <- if (median)
    inference_median_intrinsic(spobj$name, spobj$data, weights, maxiter, eps) else
    inference_mean_intrinsic(spobj$name, spobj$data, weights, maxiter, eps)
  if (!isTRUE(fit$converged))
    warning("Intrinsic location iteration did not converge (", fit$termination,
            "); returning a conditional radial estimate at the last location.",
            call. = FALSE)
  as.vector(fit[[if (median) "median" else "mean"]])
}
