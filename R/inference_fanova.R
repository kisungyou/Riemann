#' Frechet Analysis of Variance
#'
#' Compares population Frechet means and variances using the statistic of Dubey
#' and Muller (2019), equation (14). The asymptotic chi-squared calibration
#' targets equality of those summaries, not arbitrary equality of distributions.
#' The permutation version requires exchangeable observations under the stronger
#' null of identical group distributions, and targets alternatives detected by
#' the same mean-and-variance statistic.
#'
#' @param ... At least two compatible \code{riemdata} objects, each containing
#'   at least three observations. Observations must be independent within and
#'   across groups for the documented calibration; paired, repeated or clustered
#'   observations require a separate restricted resampling procedure.
#' @param maxiter Positive integer iteration budget for each Frechet mean.
#' @param eps Finite positive tolerance for the shared mean solver.
#' @param nperm Positive integer number of random label permutations.
#' @param geometry A geometry name or saved specification with a compatible
#'   mean and distance implementation. By default all groups must resolve to
#'   the same geometry.
#'
#' @details The asymptotic result requires unique population/sample means,
#'   positive variances of the squared distances to each group mean, suitable
#'   moment/metric-entropy conditions, and group proportions bounded away from
#'   zero. The cited paper gives sufficient bounded-space conditions; accepting
#'   an input manifold does not establish them for every population distribution.
#'   The function rejects failed mean fits and a degenerate empirical denominator.
#'
#'   For group proportions \eqn{\lambda_j=n_j/n}, the two denominators are
#'   \eqn{\sum_j\lambda_j/\widehat\sigma_j^2} and
#'   \eqn{\sum_j\lambda_j^2\widehat\sigma_j^2}. Distances are jointly rescaled
#'   before forming the statistic; this leaves it unchanged and avoids unit-scale
#'   overflow. Numerical rescaling does not transform the underlying geometry.
#'
#'   \code{riem.fanovaP} refits all group means for every sampled label assignment.
#'   It uses \eqn{(1+\#\{T_b\geq T_{obs}\})/(nperm+1)}, including ties. The pooled
#'   fit is unchanged by label permutation and is reused. Permutations follow R's
#'   current random-number state. The returned Monte Carlo standard error is a
#'   plug-in diagnostic for simulation variability, not inferential uncertainty
#'   in the scientific effect.
#'
#' @return An \code{htest} object with statistic, p-value, null/alternative,
#'   resolved geometry, calibration, group sizes and mean diagnostics.
#'   The permutation version additionally retains permuted statistics,
#'   \code{nperm}, and \code{mc_se}.
#'
#' @examples
#' X <- wrap.euclidean(matrix(c(-1, 0, 2, 3), ncol = 1))
#' Y <- wrap.euclidean(matrix(c(0, 1, 4, 8), ncol = 1))
#' riem.fanova(X, Y)
#' set.seed(17)
#' riem.fanovaP(X, Y, nperm = 19)
#'
#' @references
#' \insertRef{dubey_frechet_2019}{Riemann}
#' @name riem.fanova
#' @concept inference
NULL

#' @rdname riem.fanova
#' @export
riem.fanova <- function(..., maxiter = 50, eps = 1e-5, geometry = NULL) {
  prepared <- riem_fanova_prepare(list(...), geometry, maxiter, eps)
  fitted <- riem_fanova_fits(prepared)
  output <- common_fanova(fitted$pooled$distvec,
                          lapply(fitted$groups, `[[`, "distvec"),
                          prepared$geometry$manifold_id, "Supplied groups")
  output$geometry <- prepared$geometry
  output$calibration <- "asymptotic_chisquared"
  output$null.value <- c(equal_frechet_means_and_variances = 0)
  output$mean_diagnostics <- riem_fanova_diagnostics(fitted)
  output$call <- match.call()
  output
}

#' @rdname riem.fanova
#' @export
riem.fanovaP <- function(..., maxiter = 50, eps = 1e-5, nperm = 99, geometry = NULL) {
  nperm <- riem_regression_integer(nperm, "nperm", 1L)
  prepared <- riem_fanova_prepare(list(...), geometry, maxiter, eps)
  fitted <- riem_fanova_fits(prepared)
  output <- common_fanova(fitted$pooled$distvec,
                          lapply(fitted$groups, `[[`, "distvec"),
                          prepared$geometry$manifold_id, "Supplied groups")
  statistics <- numeric(nperm)
  sizes <- vapply(prepared$groups, function(x) length(x$data), integer(1))
  ends <- cumsum(sizes)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  for (b in seq_len(nperm)) {
    permutation <- sample.int(sum(sizes))
    distances <- lapply(seq_along(sizes), function(j) {
      group <- prepared$pooled
      group$data <- prepared$pooled$data[permutation[seq.int(starts[j], ends[j])]]
      fit <- riem.mean(group, geometry = prepared$geometry,
                       maxiter = prepared$maxiter, eps = prepared$eps)
      if (!isTRUE(fit$converged)) stop("A permuted group mean did not converge; no calibrated p-value is returned.",
                                      call. = FALSE)
      as.vector(basic_pdist2(group$name, group$data, list(fit$mean),
                             prepared$geometry$backend))
    })
    statistics[b] <- common_fanova(fitted$pooled$distvec, distances,
                                   prepared$geometry$manifold_id, "Permutation")$statistic
  }
  exceed <- sum(statistics >= output$statistic)
  output$p.value <- (exceed + 1) / (nperm + 1)
  output$geometry <- prepared$geometry
  output$calibration <- "random_label_permutation"
  output$null.value <- c(exchangeable_group_distributions = 0)
  output$method <- paste(output$method, "(permutation calibration)")
  output$nperm <- nperm
  output$permutation_statistics <- statistics
  output$mc_se <- sqrt(output$p.value * (1 - output$p.value) / (nperm + 1))
  output$mean_diagnostics <- riem_fanova_diagnostics(fitted)
  output$call <- match.call()
  output
}

riem_fanova_prepare <- function(groups, geometry, maxiter, eps) {
  if (length(groups) < 2L) stop("Frechet ANOVA requires at least two groups.", call. = FALSE)
  maxiter <- riem_regression_integer(maxiter, "maxiter", 1L)
  eps <- riem_legacy_positive(eps, "eps")
  specs <- lapply(groups, riem_resolve_geometry, geometry = geometry, capability = "mean")
  for (i in seq_along(groups)) {
    riem_check_newdata(groups[[1L]], groups[[i]])
    if (length(groups[[i]]$data) < 3L) stop("Each group needs at least three observations and a nondegenerate squared-distance variance.",
                                           call. = FALSE)
    if (!identical(specs[[i]], specs[[1L]])) stop("All groups must use the same geometry.", call. = FALSE)
  }
  pooled <- groups[[1L]]
  pooled$data <- unlist(lapply(groups, `[[`, "data"), recursive = FALSE)
  list(groups = groups, pooled = pooled, geometry = specs[[1L]], maxiter = maxiter, eps = eps)
}

riem_fanova_fits <- function(prepared) {
  fit <- function(data) {
    result <- riem.mean(data, geometry = prepared$geometry,
                        maxiter = prepared$maxiter, eps = prepared$eps)
    if (!isTRUE(result$converged)) stop("A Frechet mean did not converge; no calibrated test is returned.", call. = FALSE)
    result$distvec <- as.vector(basic_pdist2(data$name, data$data,
                                            list(result$mean), prepared$geometry$backend))
    result
  }
  list(pooled = fit(prepared$pooled), groups = lapply(prepared$groups, fit))
}

riem_fanova_diagnostics <- function(fits) {
  lapply(c(list(pooled = fits$pooled), fits$groups), function(fit)
    fit[c("converged", "termination", "iterations", "objective", "gradient_norm")])
}

common_fanova <- function(distall, distvecs, manifold, dataname) {
  k <- length(distvecs)
  sizes <- vapply(distvecs, length, integer(1))
  n <- length(distall)
  if (k < 2L || any(sizes < 3L) || sum(sizes) != n) {
    stop("Invalid group sizes for the Frechet ANOVA statistic.", call. = FALSE)
  }
  all_distances <- c(distall, unlist(distvecs, use.names = FALSE))
  if (any(!is.finite(all_distances)) || any(all_distances < 0)) {
    stop("Frechet ANOVA requires finite nonnegative distances.", call. = FALSE)
  }
  scale <- max(all_distances)
  if (scale == 0) stop("Squared-distance variances are degenerate; the ANOVA statistic is undefined.", call. = FALSE)
  squared <- lapply(distvecs, function(d) (d / scale)^2)
  variances <- vapply(squared, mean, numeric(1))
  sigma2 <- vapply(squared, function(d) mean((d - mean(d))^2), numeric(1))
  if (any(sigma2 <= 0)) stop("Squared-distance variances are degenerate; the ANOVA statistic is undefined.", call. = FALSE)
  lambda <- sizes / n
  pooled <- mean((distall / scale)^2)
  F <- pooled - sum(lambda * variances)
  tolerance <- 100 * .Machine$double.eps * max(pooled, variances)
  if (F < -tolerance) stop("The pooled mean objective is inconsistent with the group objectives; review the mean fits.", call. = FALSE)
  F <- max(0, F)
  U <- 0
  for (j in seq_len(k - 1L)) for (l in seq.int(j + 1L, k)) {
    U <- U + lambda[j] * lambda[l] * (variances[j] - variances[l])^2 / (sigma2[j] * sigma2[l])
  }
  terms <- c(variance = n * U / sum(lambda / sigma2),
             mean = n * F^2 / sum(lambda^2 * sigma2))
  statistic <- sum(terms)
  if (!is.finite(statistic)) stop("The ANOVA statistic is not numerically representable.", call. = FALSE)
  structure(list(statistic = c(Tn = statistic), parameter = c(df = k - 1L),
    p.value = stats::pchisq(statistic, df = k - 1L, lower.tail = FALSE),
    alternative = "at least one population Frechet mean or variance differs",
    method = paste("Frechet Analysis of Variance on", manifold), data.name = dataname,
    terms = terms, group_sizes = sizes, group_proportions = lambda,
    numerical_distance_scale = scale), class = "htest")
}
