#' Angular Central Gaussian Distribution
#' 
#' For a hypersphere \eqn{\mathcal{S}^{p-1}} in \eqn{\mathbf{R}^p}, Angular 
#' Central Gaussian (ACG) distribution \eqn{ACG_p (A)} is defined via a density
#' \deqn{f(x\vert A) = |A|^{-1/2} (x^\top A^{-1} x)^{-p/2}} 
#' with respect to the uniform measure on \eqn{\mathcal{S}^{p-1}} and \eqn{A} is 
#' a symmetric positive-definite matrix. Since \eqn{f(x\vert A) = f(-x\vert A)}, 
#' it can also be used as an axial distribution on real projective space, which is
#' unit sphere modulo \eqn{\lbrace{+1,-1\rbrace}}. One constraint we follow is that 
#' \eqn{f(x\vert A) = f(x\vert cA)} for \eqn{c > 0} in that we use a normalized 
#' version for numerical stability by restricting \eqn{tr(A)=p}.
#' 
#' @param datalist a list of length-\eqn{p} unit-norm vectors. 
#' @param A a \eqn{(p\times p)} symmetric positive-definite matrix.
#' @param n a positive integer number of samples.
#' @param log logical; return log densities when \code{TRUE}.
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
#' 
#' @details The reference measure is the uniform probability measure (total mass
#'   one), so the density equals one for the identity parameter. Density evaluation
#'   uses Cholesky factors and log determinants; \code{log=TRUE} avoids ordinary
#'   density overflow. The shape is normalized to trace \eqn{p}.
#'
#'   The fixed-point estimate carries convergence, iteration, step, and negative
#'   mean-log-likelihood attributes. Failed convergence warns. Existence and
#'   uniqueness require sufficient dispersion across all proper linear subspaces;
#'   full observed span alone does not establish these statistical conditions.
#'   Degenerate updates error, and the implementation does not claim a global
#'   existence or uniqueness certificate.
#'
#' @return 
#' \code{dacg} gives a vector of evaluated densities given samples. \code{racg} generates 
#' unit-norm vectors in \eqn{\mathbf{R}^p} wrapped in a list. \code{mle.acg} estimates 
#' the SPD matrix \eqn{A}.
#' 
#' @examples 
#' # -------------------------------------------------------------------
#' #          Example with Angular Central Gaussian Distribution
#' #
#' # Given a fixed A, generate samples and estimate A via ML.
#' # -------------------------------------------------------------------
#' ## GENERATE AND MLE in R^5
#' #  Generate data
#' Atrue = diag(5)          # true SPD matrix
#' sam1  = racg(50,  Atrue) # random samples
#' sam2  = racg(100, Atrue)
#' 
#' #  MLE
#' Amle1 = mle.acg(sam1)
#' Amle2 = mle.acg(sam2)
#' 
#' #  Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(Atrue[,5:1], axes=FALSE, main="true SPD")
#' image(Amle1[,5:1], axes=FALSE, main="MLE with n=50")
#' image(Amle2[,5:1], axes=FALSE, main="MLE with n=100")
#' par(opar)
#' 
#' @references 
#' \insertRef{tyler_statistical_1987}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @name acg
#' @concept distribution
#' @rdname acg
NULL

#' @rdname acg
#' @export
dacg <- function(datalist, A, log = FALSE) {
  data <- wrap.sphere(datalist)$data
  parameter <- riem_angular_parameter(A, "A")
  riem_angular_density(data, parameter, log)
}
#' @keywords internal
#' @noRd
dacg_internal <- function(data, A){
  n    = length(data)
  p    = length(as.vector(data[[1]]))
  coef = base::det(A)^(-0.5)
  Ainv = base::solve(A)
  
  output = rep(0,n)
  for (i in 1:n){
    tgt = as.vector(data[[i]])
    output[i] = (sum(as.vector(Ainv%*%tgt)*tgt)^(-p/2))*coef
  }
  return(output)
}

#' @rdname acg
#' @export
racg <- function(n, A) {
  n <- riem_regression_integer(n, "n", 1L)
  parameter <- riem_angular_parameter(A, "A")
  p <- nrow(parameter)
  samples <- matrix(stats::rnorm(n * p), nrow = n) %*% chol(parameter)
  lapply(seq_len(n), function(i) {
    x <- samples[i, ]
    scale <- max(abs(x))
    if (!is.finite(scale) || scale == 0) stop("The Gaussian draw cannot be normalized.", call. = FALSE)
    x <- x / scale
    x / sqrt(sum(x^2))
  })
}

#' @rdname acg
#' @export
mle.acg <- function(datalist, ...) {
  riem_angular_mle(wrap.sphere(datalist)$data, list(...), matrix_variant = FALSE)
}
