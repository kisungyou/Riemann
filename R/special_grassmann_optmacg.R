#' Estimation of Distribution Algorithm with MACG Distribution
#' 
#' For a function \eqn{f : Gr(k,p) \rightarrow \mathbf{R}}, find the minimizer 
#' and the attained minimum value with estimation of distribution algorithm 
#' using MACG distribution.
#' 
#' @param func a function to be \emph{minimized}.
#' @param p dimension parameter as in \eqn{Gr(k,p)}.
#' @param k dimension parameter as in \eqn{Gr(k,p)}.
#' @param ... extra parameters including\describe{
#' \item{n.start}{number of runs; algorithm is executed \code{n.start} times (default: 10).}
#' \item{maxiter}{maximum number of iterations for each run (default: 100).}
#' \item{popsize}{the number of samples generated at each step for stochastic search (default: 100).}
#' \item{ratio}{ratio in \eqn{(0,1)} where top \code{ratio*popsize} samples are chosen for parameter update (default: 0.25).}
#' \item{print.progress}{a logical; if \code{TRUE}, it prints each iteration (default: \code{FALSE}).}
#' }
#' 
#' @return a named list containing: \describe{
#' \item{cost}{minimized function value.}
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
#' A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
#' myfunc <- function(p){                # cost function to minimize
#'   return(sum(-diag(t(p)%*%A%*%p)))
#' } 
#' 
#' ## SOLVE THE OPTIMIZATION PROBLEM
#' Aout = grassmann.optmacg(myfunc, p=5, k=3, popsize=100, n.start=30)
#' 
#' ## COMPUTE EIGENVALUES
#' #  1. USE SOLUTIONS TO THE ABOVE OPTIMIZATION 
#' abase   = Aout$solution
#' eig3sol = sort(diag(t(abase)%*%A%*%abase), decreasing=TRUE)
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
grassmann.optmacg <- function(func, p, k, ...){
  # Preprocessing
  # 1. function
  FNAME = "grassmann.optmacg"
  if (!is.function(func)){
    stop(paste0("* ",FNAME," : an input 'func' should be a function."))
  }
  # 2. implicit parameters : maxiter, popsize, ratio
  params = list(...)
  pnames = names(params)
  
  nstart  = ifelse(("n.start"%in%pnames), max(10, params$n.start), 10) 
  maxiter = ifelse(("maxiter"%in%pnames), max(20, round(params$maxiter)), 100)
  popsize = ifelse(("popsize"%in%pnames), max(50, round(params$popsize)), 100); popsize = max(popsize, round(p*k)*3)
  ratio   = ifelse(("ratio"%in%pnames), params$ratio, 0.25)
  printer = ifelse(("print.progress"%in%pnames), as.logical(params$print.progress), FALSE)
  
  if ((ratio <= 0)||(ratio >=1)){
    stop(paste0("* ",FNAME," : 'ratio' should be a number in (0,1)."))
  }
  myp = round(p)
  myk = round(k)
  toppick = round(ratio*popsize)
  
  ## MAIN COMPUTATION
  #  set up initial covariance matrices
  start_Sigma = list()
  start_Sigma[[1]] = base::diag(myp)
  for (j in 2:nstart){
    sam_mats  = runif_stiefel(myp,myk,2*popsize) # 3d array with C++
    vec_fvals = rep(0,(2*popsize))
    for (i in 1:(2*popsize)){
      vec_fvals[i] = as.double(func(as.matrix(sam_mats[,,i])))
    }
    min_ids = base::order(vec_fvals)[1:toppick]
    min_mat = list()
    for (i in 1:toppick){
      min_mat[[i]] = as.matrix(sam_mats[,,min_ids[i]])
    }
    tmpsig = mle.macg(min_mat)
    start_Sigma[[j]] = as.matrix(Matrix::nearPD((tmpsig+t(tmpsig))/2)$mat)
  }
  #  compute
  single_runs = list()
  single_vals = rep(0,nstart)
  for (i in 1:nstart){
    single_runs[[i]] = grassmann.optmacg.single(func, myp, myk, start_Sigma[[i]],
                                                maxiter, popsize, toppick, i, printer)
    single_vals[i] = single_runs[[i]]$cost
  }
  

  ## RETURN
  opt_id  = which.min(single_vals)
  opt_run = single_runs[[opt_id]]
  return(opt_run)
}
#' @keywords internal
#' @noRd
grassmann.optmacg.single <- function(func, myp, myk, Sigma, maxiter, popsize, toppick, runid, printer){
  mleS    = Sigma
  min_f   = 10000000
  min_mat = NULL
  
  counter = 0
  for (it in 1:maxiter){
    # step 1. random generation
    random_sample = rmacg(popsize, myk, mleS)
    # step 2. evaluate function values
    vec_fvals = rep(0,popsize)
    for (i in 1:popsize){
      vec_fvals[i] = as.double(func(random_sample[[i]]))
    }
    # step 3. update the minimal ones
    if (min(vec_fvals) < min_f){
      idmin   = which.min(vec_fvals)
      min_f   = vec_fvals[idmin]
      min_mat = random_sample[[idmin]]
      counter = 0
    } else {
      counter = counter + 1
    }
    if (counter >= 5){
      break
    }
    # step 4. update mleS
    top_ids  = base::order(vec_fvals)[1:toppick]
    top_mats = random_sample[top_ids]
    mleS     = mle.macg(top_mats)
    mleS     = as.matrix(Matrix::nearPD((mleS + t(mleS))/2)$mat)
    if (printer){
      print(paste0("* grassmann.optmacg : run ",runid," & iter ",it,"/",maxiter,"; current fmin value=",round(min_f,5),"."))
    }
  }
  output = list()
  output$cost     = min_f
  output$solution = min_mat
  return(output)
}

