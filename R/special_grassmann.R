## SOME SPECIAL FUNCTIONS ON GRASSMANN 
#  (01) grassmann.runif
#  (02) grassmann.utest
#  (03) grassmann.optmacg

#  (03) grassmann.optmacg ======================================================
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


#  (01) grassmann.runif ========================================================
#' Generate Uniform Samples on Grassmann Manifold
#' 
#' It generates \eqn{n} random samples from Grassmann manifold \eqn{Gr(k,p)}.
#' 
#' @param n number of samples to be generated.
#' @param k dimension of the subspace.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of \eqn{(p\times k)} basis of \eqn{k}-subspaces.}
#' \item{\code{"array"}}{a \eqn{(p\times k\times n)} 3D array whose slices are \eqn{k}-subspace basis.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.grassmann}} for more details.}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' @seealso \code{\link{stiefel.runif}}, \code{\link{wrap.grassmann}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Draw Samples on Grassmann Manifold 
#' #-------------------------------------------------------------------
#' #  Multiple Return Types with 3 Observations of 5-dim subspaces in R^10
#' dat.list = grassmann.runif(n=3, k=5, p=10, type="list")
#' dat.arr3 = grassmann.runif(n=3, k=5, p=10, type="array")
#' dat.riem = grassmann.runif(n=3, k=5, p=10, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @concept grassmann
#' @export
grassmann.runif <- function(n, k, p, type=c("list","array","riemdata")){
  #  PREPROCESSING
  n = round(n)
  p = round(p) # k-frame in R^p
  k = round(k)
  if (k > p){
    stop("* grassmann.runif : 'k <= p' is a required condition.")
  }
  retype = ifelse(missing(type),"riemdata",
                  match.arg(tolower(type), c("list","array","riemdata")))
  
  #  GENERATE, WRAP, RETURN
  tmpout = runif_stiefel(p,k,n) # C++ version
  if (all(retype=="array")){
    return(tmpout)
  } else {
    tmpobj = wrap.grassmann(tmpout)
    if (all(retype=="list")){
      return(tmpobj$data)
    } else {
      return(tmpobj)
    }
  }
}

#  (02) grassmann.utest ========================================================
#' Test of Uniformity on Grassmann Manifold
#' 
#' Given the data on Grassmann manifold \eqn{Gr(k,p)}, it tests whether the 
#' data is distributed uniformly.
#' 
#' @param grobj a S3 \code{"riemdata"} class of Grassmann-valued data.
#' @param method (case-insensitive) name of the test method containing \describe{
#' \item{\code{"Bing"}}{Bingham statistic.}
#' \item{\code{"BingM"}}{modified Bingham statistic with better order of error.}
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
#' @seealso \code{\link{wrap.grassmann}}
#' @examples 
#' #-------------------------------------------------------------------
#' #   Compare Bingham's original and modified versions of the test
#' # 
#' # Test 1. sample uniformly from Gr(2,4)
#' # Test 2. use perturbed principal components from 'iris' data in R^4
#' #         which is concentrated around a point to reject H0.
#' #-------------------------------------------------------------------
#' ## Data Generation
#' #  1. uniform data
#' myobj1 = grassmann.runif(n=100, k=2, p=4)
#' 
#' #  2. perturbed principal components
#' data(iris)
#' irdat = list()
#' for (n in 1:100){
#'    tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
#'    irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
#' }
#' myobj2 = wrap.grassmann(irdat)
#' 
#' ## Test 1 : uniform data
#' grassmann.utest(myobj1, method="Bing")
#' grassmann.utest(myobj1, method="BingM")
#' 
#' ## Tests : iris data
#' grassmann.utest(myobj2, method="bINg")   # method names are 
#' grassmann.utest(myobj2, method="BiNgM")  # CASE - INSENSITIVE !
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @concept grassmann
#' @export
grassmann.utest <- function(grobj, method=c("Bing","BingM")){
  #  CHECK INPUT
  DNAME    = deparse(substitute(grobj)) # borrowed from HDtest
  FNAME    = "grassmann.utest"
  check_inputmfd(grobj, FNAME)
  mymethod = ifelse(missing(method),"bing",
                    match.arg(tolower(method),c("bing","bingm")))
  
  #  COMPUTE AND RETURN
  output <- switch(mymethod,
                   bing  = gr.utest.bing(grobj, DNAME, is.modified = FALSE),
                   bingm = gr.utest.bing(grobj, DNAME, is.modified = TRUE))
  return(output)
}
#' @keywords internal
#' @noRd
gr.utest.bing <- function(grobj, dname, is.modified=TRUE){
  # PREPARE
  p = grobj$size[1] # parameters
  r = grobj$size[2]
  n = length(grobj$data)
  
  # COMPUTE THE MEAN
  Ybar = array(0,c(p,p))
  for (i in 1:n){
    Ytgt = grobj$data[[i]]
    Ybar = Ybar + (Ytgt%*%t(Ytgt))/n
  }
  
  # COMPUTE THE STATISTIC
  S = (((p-1)*p*(p+2))/(2*r*(p-r)))*n*(sum(diag(Ybar%*%Ybar)) - ((r^2)/p))
  
  # BRANCHING
  if (is.modified){
    p2 = p*p
    B0 = ((2*p2*(p-1)*(p+2))-(r*(p-r)*(5*p2 + 2*p + 8)))/(6*r*(p-r)*(p-2)*(p+4))
    B1 = (-(4*p2*(p-1)*(p+2))+(r*(p-r)*((13*p2)+(10*p)-8)))/(3*r*(p-r)*(p2+p+2)*(p-2)*(p+4))
    B2 = (4*((p-(2*r))^2)*(p2+p-2))/(3*r*(p-r)*(p-2)*(p+4)*(p2+p+2)*(p2+p+6))
    
    thestat = S*(1-(1/n)*(B0 + B1*S + B2*(S^2)))
    thedf   = round((p-1)*(p+2)/2)
    hname   = "Modified Bingham Test of Uniformity on Grassmann Manifold"
    pvalue  = stats::pchisq(thestat, df=thedf, lower.tail=FALSE)
  } else {
    hname   = "Bingham Test of Uniformity on Grassmann Manifold"
    thestat = S
    thedf   = round((p-1)*(p+2)/2)
    pvalue  = stats::pchisq(thestat, df=thedf, lower.tail=FALSE)
  }
  
  # DETERMINATION
  Ha      = paste("data is not uniformly distributed on Gr(",r,",",p,").",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}