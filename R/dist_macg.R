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
#' @param n the number of samples to be generated.
#' @param r the number of basis.
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
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
dmacg <- function(datalist, Sigma){
  ## CHECK INPUT
  FNAME = "dmacg"
  myobj = wrap.stiefel(datalist)
  mysig = as.matrix(Sigma)
  if (!check_spdmat(mysig)){
    stop(paste0("* ",FNAME," : 'Sigma' should be a symmetric positive-definite matrix."))
  }
  
  ## COMPUTE
  output = macg_density(myobj$data, mysig)
  return(as.vector(output))
}

#' @rdname macg
#' @export
rmacg <- function(n, r, Sigma){
  ## INITIALIZATION
  FNAME = "rmacg"
  mysig = as.matrix(Sigma)
  if (!check_spdmat(mysig)){
    stop(paste0("* ",FNAME," : 'Sigma' should be a symmetric positive-definite matrix."))
  }
  myp   = base::nrow(mysig)
  myn   = max(1, round(n))
  mysig = (mysig/sum(diag(mysig)))*myp
  myr   = max(round(r), 1)
  if (myr > myp){
    stop(paste0("* ",FNAME," : we require 'r<=nrow(Sigma)'."))
  }

  ## COMPUTE
  outcube = macg_sample(myn, myr, mysig)
  output  = list()
  for (i in 1:n){
    output[[i]] = as.matrix(outcube[,,i])
  }
  return(output)
}

#' @rdname macg
#' @export
mle.macg <- function(datalist, ...){
  ## CHECK INPUT
  myobj = wrap.stiefel(datalist)
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 50))
  myeps  = min(1e-5, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-5)))
  
  ## COMPUTE AND RETURN
  output = macg_mle(myobj$data, myiter, myeps)
  return(output)
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