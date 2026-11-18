#' Estimation of Distribution Algorithm with MACG Distribution
#' 
#' Search for a subspace minimizing a basis-invariant objective by adapting a
#' matrix angular central Gaussian (MACG) sampling distribution to elite frames.
#' 
#' @param func a function to be \emph{minimized}, returning one finite numeric value.
#' @param p dimension parameter as in \eqn{Gr(k,p)}.
#' @param k dimension parameter as in \eqn{Gr(k,p)}.
#' @param ... named controls including\describe{
#' \item{n.start}{number of runs; algorithm is executed \code{n.start} times (default: 10).}
#' \item{maxiter}{maximum number of iterations for each run (default: 100).}
#' \item{popsize}{the number of samples generated at each step for stochastic search (default: 100).}
#' \item{ratio}{ratio in \eqn{(0,1)} where top \code{ratio*popsize} samples are chosen for parameter update (default: 0.25).}
#' \item{print.progress}{a logical; if \code{TRUE}, it prints each iteration (default: \code{FALSE}).}
#' }
#' 
#' @details The first run starts with uniform samples. Each additional run
#' begins with a uniform pilot population twice the requested size. At each
#' update, the best frames determine a trace-normalized MACG proposal via
#' at most 50 fixed-point steps. Each step mixes in a fraction
#' \code{sqrt(.Machine$double.eps)} of the identity to keep the search
#' distribution positive definite even when the elite subspaces coincide.
#' This regularization applies only to the optimizer's sampling proposal.
#' A run ends after five generations without improvement or at \code{maxiter}.
#' For \eqn{k=p}, the Grassmannian contains only one subspace.
#'
#' @return a named list containing: \describe{
#' \item{cost}{smallest function value attained during the iterative search runs; a global minimum is not guaranteed.}
#' \item{solution}{a \eqn{(p\times k)} matrix that attains the \code{cost}.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #               Optimization for Eigen-Decomposition
#' #
#' # Given (5x5) covariance matrix S, eigendecomposition is can be 
#' # considered as an optimization on Grassmann manifold. Here, 
#' # we are trying to find top 3 eigenvalues and compare.
#' #-------------------------------------------------------------------
#' \donttest{
#' ## PREPARE
#' set.seed(17)
#' A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
#' myfunc <- function(p){                # cost function to minimize
#'   return(sum(-diag(t(p)%*%A%*%p)))
#' } 
#' 
#' ## SOLVE THE OPTIMIZATION PROBLEM
#' Aout = grassmann.optmacg(myfunc, p=5, k=3, popsize=100, n.start=3, maxiter=30)
#' 
#' ## COMPUTE EIGENVALUES
#' #  1. USE SOLUTIONS TO THE ABOVE OPTIMIZATION 
#' abase   = Aout$solution
#' eig3sol = eigen(crossprod(abase, A %*% abase), symmetric=TRUE)$values
#'
#' #  2. USE BASIC 'EIGEN' FUNCTION
#' eig3dec = sort(eigen(A)$values, decreasing=TRUE)[1:3]
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' yran = c(min(min(eig3sol),min(eig3dec))*0.95,
#'          max(max(eig3sol),max(eig3dec))*1.05)
#' plot(1:3, eig3sol, type="b", col="red",  pch=19, ylim=yran,
#'      xlab="index", ylab="eigenvalue", main="compare top 3 eigenvalues")
#' lines(1:3, eig3dec, type="b", col="blue", pch=19)
#' legend(1.55, max(yran), legend=c("optimization","decomposition"), col=c("red","blue"),
#'        lty=rep(1,2), pch=19)
#' par(opar)
#' }
#' 
#' @concept grassmann
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
grassmann.optmacg <- function(func, p, k, ...) {
  if (!is.function(func)) stop("'func' must be a function.", call. = FALSE)
  p <- riem_regression_integer(p, "p", 1L)
  k <- riem_regression_integer(k, "k", 1L, p)
  controls <- list(n.start = 10L, maxiter = 100L, popsize = 100L, ratio = 0.25,
                   print.progress = FALSE)
  supplied <- riem_legacy_parameters(list(...), names(controls))
  controls[names(supplied)] <- supplied
  nstart <- riem_regression_integer(controls$n.start, "n.start", 1L)
  maxiter <- riem_regression_integer(controls$maxiter, "maxiter", 1L)
  popsize <- riem_regression_integer(controls$popsize, "popsize", 1L,
                                     .Machine$integer.max %/% 2L)
  ratio <- riem_legacy_positive(controls$ratio, "ratio")
  if (ratio >= 1) stop("'ratio' must be in (0, 1).", call. = FALSE)
  if (!is.logical(controls$print.progress) || length(controls$print.progress) != 1L ||
      is.na(controls$print.progress)) {
    stop("'print.progress' must be TRUE or FALSE.", call. = FALSE)
  }
  elite_count <- max(1L, round(ratio * popsize))
  if (k == p) {
    frame <- diag(p)
    return(list(cost = grassmann_optmacg_value(func, frame), solution = frame))
  }

  best <- list(cost = Inf, solution = NULL)
  for (run in seq_len(nstart)) {
    proposal <- diag(p)
    run_best <- Inf
    stale <- 0L
    # Generation zero initializes additional starts from a larger uniform pilot.
    for (generation in seq.int(as.integer(run == 1L), maxiter)) {
      size <- if (generation == 0L) 2L * popsize else popsize
      frames <- rmacg(size, k, proposal)
      costs <- vapply(frames, function(x) grassmann_optmacg_value(func, x), numeric(1))
      ranking <- order(costs)
      winner <- ranking[1L]
      if (costs[winner] < best$cost) {
        best <- list(cost = unname(costs[winner]), solution = frames[[winner]])
      }
      if (costs[winner] < run_best) {
        run_best <- costs[winner]
        stale <- 0L
      } else {
        stale <- stale + 1L
      }
      if (controls$print.progress && generation > 0L) {
        message(sprintf("grassmann.optmacg: run %d, iteration %d/%d, best cost %.6g",
                        run, generation, maxiter, run_best))
      }
      if (stale >= 5L || generation == maxiter) break
      proposal <- grassmann_optmacg_proposal(frames[ranking[seq_len(elite_count)]], proposal)
    }
  }
  best
}

# Fit a proposal, not an unregularized MLE: a single or repeatedly selected
# elite subspace can legitimately have zero variance in some directions.
grassmann_optmacg_proposal <- function(frames, proposal) {
  p <- nrow(proposal)
  ridge <- sqrt(.Machine$double.eps)
  for (i in seq_len(50L)) {
    update <- riem_angular_update(frames, proposal)
    update <- update * ((1 - ridge) * p / sum(diag(update)))
    diag(update) <- diag(update) + ridge
    step <- norm(update - proposal, "F")
    proposal <- update
    if (step <= 1e-5) break
  }
  proposal
}

#' @keywords internal
#' @noRd
grassmann_optmacg_value <- function(func, x) {
  value <- func(x)
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value))
    stop("* grassmann.optmacg : 'func' must return a finite numeric scalar.", call. = FALSE)
  as.double(value)
}
