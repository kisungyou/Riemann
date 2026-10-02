#' Two-Sample Test with Wasserstein Metric
#' 
#' Given \eqn{M} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}} and 
#' \eqn{N} observations \eqn{Y_1, Y_2, \ldots, Y_N \in \mathcal{M}}, permutation 
#' test based on the Wasserstein metric (see \code{\link{riem.wasserstein}} for 
#' more details) is applied to test whether two distributions are same or not, i.e.,
#' \deqn{H_0~:~\mathcal{P}_X = \mathcal{P}_Y}
#' with Wasserstein metric \eqn{\mathcal{W}_p} being the measure of discrepancy 
#' between two samples.
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param p an exponent for Wasserstein distance \eqn{\mathcal{W}_p} (default: 2).
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{nperm}{the number of permutations (default: 999).}
#' \item{use.smooth}{a logical, default \code{FALSE}. \code{TRUE} requests an experimental IPOT approximation, without an optimality certificate.}
#' }
#' 
#' @details Random label permutations require independent observations that are
#'   exchangeable under the identical-distributions null. This procedure does not
#'   handle paired, clustered or repeated observations. The Monte Carlo p-value
#'   includes ties and uses the plus-one correction. All assignments use the same
#'   solver and exponent. Saved geometry, permuted statistics, and the Monte Carlo
#'   standard-error diagnostic are retained.
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
#' #  it is expected to return a very small number, but 
#' #  small number of 'nperm' may not give a reasonable p-value.
#' \donttest{
#' riem.test2wass(myriem1, myriem2, nperm=99, use.smooth=FALSE)
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
#'   pvals[i] = riem.test2wass(Xriem, Yriem, nperm=999)$p.value
#'   print(paste0("iteration ",i,"/",ntest," complete.."))
#' }
#' 
#' emperr = round(sum((pvals <= 0.05))/ntest, 5)
#' print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
#' }
#' 
#' @concept inference
#' @export
riem.test2wass <- function(riemobj1, riemobj2, p = 2, geometry = NULL, ...) {
  inputs <- riem_transport_inputs(riemobj1, riemobj2, p, geometry)
  parameters <- riem_legacy_parameters(list(...), c("nperm", "use.smooth"))
  nperm <- if (is.null(parameters$nperm)) 999L else
    riem_regression_integer(parameters$nperm, "nperm", 1L)
  smooth <- if (is.null(parameters$use.smooth)) FALSE else parameters$use.smooth
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    stop("use.smooth must be TRUE or FALSE.", call. = FALSE)
  }
  if (smooth) warning("use.smooth uses an experimental IPOT approximation without an optimality certificate; the statistic is not an exact Wasserstein distance.",
                       call. = FALSE)
  M <- length(riemobj1$data)
  N <- length(riemobj2$data)
  wx <- rep(1 / M, M)
  wy <- rep(1 / N, N)
  distances <- basic_pdist(riemobj1$name, c(riemobj1$data, riemobj2$data), inputs$geometry$backend)
  compute <- function(ix, iy) riem_transport_solve(distances[ix, iy, drop = FALSE],
                                                  inputs$p, wx, wy, smooth)$distance
  observed <- compute(seq_len(M), M + seq_len(N))
  permuted <- numeric(nperm)
  for (b in seq_len(nperm)) {
    ids <- sample.int(M + N)
    permuted[b] <- compute(ids[seq_len(M)], ids[M + seq_len(N)])
  }
  pvalue <- (1 + sum(permuted >= observed)) / (nperm + 1)
  structure(list(statistic = c(Wmn = observed), p.value = pvalue,
    alternative = "the group distributions differ",
    null.value = c(exchangeable_group_distributions = 0),
    method = if (smooth) "Two-sample permutation test using approximate transport costs" else
      "Wasserstein two-sample permutation test",
    data.name = "Supplied samples", geometry = inputs$geometry, p = inputs$p,
    calibration = "random_label_permutation", nperm = nperm,
    permutation_statistics = permuted,
    mc_se = sqrt(pvalue * (1 - pvalue) / (nperm + 1)),
    solver = if (smooth) "IPOT_approximation" else "linear_program"), class = "htest")
}
