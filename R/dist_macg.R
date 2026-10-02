#' Matrix Angular Central Gaussian Distribution
#' 
#' For Stiefel and Grassmann manifolds \eqn{St(r,p)} and \eqn{Gr(r,p)}, the matrix 
#' variant of ACG distribution is known as Matrix Angular Central Gaussian (MACG) 
#' distribution \eqn{MACG_{p,r}(\Sigma)} with density
#' \deqn{f(X\vert \Sigma) = |\Sigma|^{-r/2} |X^\top \Sigma^{-1} X|^{-p/2}}
#' where \eqn{\Sigma} is a \eqn{(p\times p)} symmetric positive-definite matrix. 
#' Similar to vector-variate ACG case, we follow a convention that \eqn{tr(\Sigma)=p}.
#' 
#' @param datalist a list of \eqn{(p\times r)} orthonormal matrices.
#' @param Sigma a \eqn{(p\times p)} symmetric positive-definite matrix.
#' @param n a positive integer number of samples.
#' @param log logical; return log densities when \code{TRUE}.
#' @param r the number of basis.
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
#' 
#' @details The density is relative to the invariant uniform probability measure
#'   on the Stiefel manifold, or its push-forward on the Grassmann quotient, and
#'   equals one at \eqn{\Sigma=I}. Cholesky/log-determinant evaluation avoids
#'   determinant products that overflow at extreme parameter scales.
#'
#'   Sampling uses the polar factor of a Gaussian matrix. The fixed-point estimate
#'   is trace normalized and carries convergence, iteration, step, and negative
#'   mean-log-likelihood attributes. Failed convergence warns; degenerate updates
#'   error. Full observed span does not certify the subspace-dispersion conditions
#'   needed for an identifiable, unique interior estimate. When \eqn{r=p>1}, the
#'   density is uniform regardless of \eqn{\Sigma}, so estimation errors.
#'
#' @return 
#' \code{dmacg} gives a vector of evaluated densities given samples. \code{rmacg} generates  
#' \eqn{(p\times r)} orthonormal matrices wrapped in a list. \code{mle.macg} estimates 
#' the SPD matrix \eqn{\Sigma}.
#' 
#' @examples 
#' # -------------------------------------------------------------------
#' #          Example with Matrix Angular Central Gaussian Distribution
#' #
#' # Given a fixed Sigma, generate samples and estimate Sigma via ML.
#' # -------------------------------------------------------------------
#' ## GENERATE AND MLE in St(2,5)/Gr(2,5)
#' #  Generate data
#' Strue = diag(5)                  # true SPD matrix
#' sam1  = rmacg(n=50,  r=2, Strue) # random samples
#' sam2  = rmacg(n=100, r=2, Strue) # random samples
#' 
#' #  MLE
#' Smle1 = mle.macg(sam1)
#' Smle2 = mle.macg(sam2)
#' 
#' #  Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(Strue[,5:1], axes=FALSE, main="true SPD")
#' image(Smle1[,5:1], axes=FALSE, main="MLE with n=50")
#' image(Smle2[,5:1], axes=FALSE, main="MLE with n=100")
#' par(opar)
#' 
#' @seealso \code{\link{acg}}
#' 
#' @references 
#' \insertRef{chikuse_matrix_1990}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' Kent JT, Ganeiber AM, Mardia KV (2013). "A new method to simulate the Bingham and related distributions in directional data analysis with applications." \emph{arXiv:1310.8110}.
#' 
#' @name macg
#' @concept distribution
#' @rdname macg
NULL

#' @rdname macg
#' @export
dmacg <- function(datalist, Sigma, log = FALSE) {
  data <- wrap.stiefel(datalist)$data
  parameter <- riem_angular_parameter(Sigma, "Sigma")
  riem_angular_density(data, parameter, log)
}

#' @rdname macg
#' @export
rmacg <- function(n, r, Sigma) {
  n <- riem_regression_integer(n, "n", 1L)
  parameter <- riem_angular_parameter(Sigma, "Sigma")
  p <- nrow(parameter)
  r <- riem_regression_integer(r, "r", 1L, p)
  root <- t(chol(parameter))
  lapply(seq_len(n), function(i) {
    raw <- root %*% matrix(stats::rnorm(p * r), p, r)
    polar <- svd(raw, nu = r, nv = r)
    if (min(polar$d) <= 0) stop("The Gaussian draw has deficient column rank.", call. = FALSE)
    polar$u %*% t(polar$v)
  })
}

#' @rdname macg
#' @export
mle.macg <- function(datalist, ...) {
  riem_angular_mle(wrap.stiefel(datalist)$data, list(...), matrix_variant = TRUE)
}

# A = matrix(runif(100*5),ncol=5)
# S = t(A)%*%A
# S = diag(5)
# 
# sam1 = rmacg(100, 3, S); mle1 = mle.macg(sam1)
# sam2 = rmacg(100, 3, S); mle2 = mle.macg(sam2)
# 
# par(mfrow=c(1,3),pty="s")
# image(S, axes=FALSE, main="true")
# image(mle1, axes=FALSE, main="MLE 1")
# image(mle2, axes=FALSE, main="MLE 2")