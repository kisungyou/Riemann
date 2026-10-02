#' Two-Sample Test modified from Biswas and Ghosh (2014)
#' 
#' Given \eqn{M} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}} and 
#' \eqn{N} observations \eqn{Y_1, Y_2, \ldots, Y_N \in \mathcal{M}}, perform the permutation test of equal distribution
#' \deqn{H_0~:~\mathcal{P}_X = \mathcal{P}_Y}
#' by the method from Biswas and Ghosh (2014). The method, originally proposed 
#' for Euclidean-valued data, is adapted to the general Riemannian manifold 
#' with intrinsic/extrinsic distance. 
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry A geometry name or saved specification.
#' @details The statistic compares within-group and cross-group mean distances.
#'   Random-label calibration requires exchangeability under the equal-distribution
#'   null, with independent observations; it is not valid for unaccounted paired
#'   or repeated observations. Generalizing the distance to a manifold does not
#'   establish consistency against every alternative. Monte Carlo p-values include
#'   the observed assignment and all ties, and cannot be zero. Each group must
#'   contain at least two observations.
#' @param ... extra parameters including\describe{
#' \item{nperm}{the number of permutations (default: 999).}
#' }
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with two types
#' #
#' # class 1 : 20 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 30 perturbed data points near (0,1,0) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata1 = list()
#' mydata2 = list()
#' for (i in 1:20){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata1[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 1:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem1 = wrap.sphere(mydata1)
#' myriem2 = wrap.sphere(mydata2)
#' 
#' ## PERFORM PERMUTATION TEST
#' #  it is expected to return a very small number.
#' \donttest{
#' riem.test2bg14(myriem1, myriem2, nperm=999)
#' }
#' 
#' \dontrun{
#' ## CHECK WITH EMPIRICAL TYPE-1 ERROR
#' set.seed(777)
#' ntest = 1000
#' pvals = rep(0,ntest)
#' 
#' for (i in 1:ntest){
#'   X = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Y = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Xnorm = X/sqrt(rowSums(X^2))
#'   Ynorm = Y/sqrt(rowSums(Y^2))
#'   
#'   Xriem = wrap.sphere(Xnorm)
#'   Yriem = wrap.sphere(Ynorm)
#'   pvals[i] = riem.test2bg14(Xriem, Yriem, nperm=999)$p.value
#' }
#' 
#' emperr = round(sum((pvals <= 0.05))/ntest, 5)
#' print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
#' }
#' 
#' @references
#' \insertRef{biswas_nonparametric_2014a}{Riemann}
#' 
#' \insertRef{you_revisiting_2020a}{Riemann}
#' 
#' @concept inference
#' @export
riem.test2bg14 <- function(riemobj1, riemobj2, geometry = NULL, ...) {
  riem_validate_data(riemobj1)
  riem_check_newdata(riemobj1, riemobj2)
  spec <- riem_resolve_geometry(riemobj1, geometry, capability = "distance")
  other <- riem_resolve_geometry(riemobj2, geometry, capability = "distance")
  if (!identical(spec, other)) stop("Both groups must use the same geometry.", call. = FALSE)
  m <- length(riemobj1$data)
  n <- length(riemobj2$data)
  if (m < 2L || n < 2L) stop("Each group requires at least two observations.", call. = FALSE)
  options <- list(...)
  if (length(options) && (is.null(names(options)) || any(!names(options) %in% "nperm") ||
                          anyDuplicated(names(options)))) stop("Only the named nperm option is supported.", call. = FALSE)
  nperm <- riem_regression_integer(if (is.null(options$nperm)) 999L else options$nperm,
                                    "nperm", 1L)
  distances <- basic_pdist(riemobj1$name, c(riemobj1$data, riemobj2$data), spec$backend)
  if (any(!is.finite(distances)) || any(distances < 0)) stop("Invalid pairwise distances.", call. = FALSE)
  statistic <- function(ix) {
    iy <- setdiff(seq_len(m + n), ix)
    R_eqdist_2014BG_statistic(distances[ix, ix, drop = FALSE], distances[iy, iy, drop = FALSE],
                             distances[ix, iy, drop = FALSE])
  }
  observed <- statistic(seq_len(m))
  permuted <- replicate(nperm, statistic(sample.int(m + n, m)))
  pvalue <- (1 + sum(permuted >= observed)) / (nperm + 1)
  structure(list(statistic = c(Tmn = observed), p.value = pvalue,
    alternative = "the group distributions differ in the distance summaries",
    null.value = c(exchangeable_group_distributions = 0),
    method = "Biswas-Ghosh distance-summary random-label test",
    data.name = paste(deparse(substitute(riemobj1)), "and", deparse(substitute(riemobj2))),
    geometry = spec, calibration = "random_label_permutation", nperm = nperm,
    permutation_statistics = permuted,
    mc_se = sqrt(pvalue * (1 - pvalue) / (nperm + 1))), class = "htest")
}

R_eqdist_2014BG_statistic <- function(DX, DY, DXY) {
  m <- nrow(DXY)
  n <- ncol(DXY)
  if (m < 2L || n < 2L) stop("Distance-summary statistic requires two observations per group.", call. = FALSE)
  within_x <- mean(DX[upper.tri(DX)])
  within_y <- mean(DY[upper.tri(DY)])
  between <- mean(DXY)
  value <- (within_x - between)^2 + (within_y - between)^2
  if (!is.finite(value)) stop("The distance-summary statistic is nonfinite; review units.", call. = FALSE)
  value
}
