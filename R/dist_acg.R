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
#' @param n the number of samples to be generated.
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
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
dacg <- function(datalist, A){
  ## INITIALIZATION
  FNAME = "dacg"
  myobj = wrap.sphere(datalist)
  myA   = as.matrix(A)
  if (!check_spdmat(myA)){
    stop(paste0("* ",FNAME," : 'A' should be a symmetric positive-definite matrix."))
  }
  if (base::nrow(myA)!=length(as.vector(myobj$data[[1]]))){
    stop(paste0("* ",FNAME," : 'A' should have ",length(as.vector(myobj$data[[1]]))," columns and rows."))
  }
  
  ## COMPUTATION
  output = acg_density(myobj$data, myA)
  # output = dacg_internal(myobj$data, myA)
  
  # if (TRUE){ # adjust with surface measure
  #   p = base::nrow(myA)
  #   Cp = (2*(pi^(p/2)))/base::gamma(p/2)
  #   output = as.vector(output)/Cp
  # }
  return(as.vector(output))
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
racg <- function(n, A){
  ## INITIALIZATION
  FNAME = "racg"
  myA   = as.matrix(A)
  if (!check_spdmat(myA)){
    stop(paste0("* ",FNAME," : 'A' should be a symmetric positive-definite matrix."))
  }
  myp = base::nrow(myA)
  myn = max(1, round(n))
  myA = (myA/sum(diag(myA)))*myp # normalize to "tr(A)=p"
  
  ## COMPUTE
  cppsam = cpp_rmvnorm(myn, rep(0,myp), myA)
  output = list()
  for (i in 1:myn){
    tgt         = as.vector(cppsam[i,])
    output[[i]] = tgt/sqrt(sum(tgt^2))
  }
  return(output)
}

#' @rdname acg
#' @export
mle.acg <- function(datalist, ...){
  ## INITIALIZATION
  myobj  = wrap.sphere(datalist)
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 50))
  myeps  = min(1e-5, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-5)))
  
  ## COMPUTE AND RETURN
  output = acg_mle(myobj$data, myiter, myeps)
  return(output)
}