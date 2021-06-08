## SOME SPECIAL FUNCTIONS ON STIEFEL 
#  (01) stiefel.runif
#  (02) stiefel.utest
#  (03) stiefel.optSA

#  (01) stiefel.runif ==========================================================
#' Generate Uniform Samples on Stiefel Manifold
#' 
#' It generates \eqn{n} random samples from Stiefel manifold \eqn{St(k,p)}.
#' 
#' @param n number of samples to be generated.
#' @param k dimension of the frame.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of \eqn{(p\times k)} basis of \eqn{k}-frames.}
#' \item{\code{"array"}}{a \eqn{(p\times k\times n)} 3D array whose slices are \eqn{k}-frame basis.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.stiefel}} for more details.}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' @seealso \code{\link{wrap.stiefel}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Draw Samples on Stiefel Manifold 
#' #
#' # Try Different Return Types with 3 Observations of 5-frames in R^10
#' #-------------------------------------------------------------------
#' #  GENERATION
#' dat.list = stiefel.runif(n=3, k=5, p=10, type="list")
#' dat.arr3 = stiefel.runif(n=3, k=5, p=10, type="array")
#' dat.riem = stiefel.runif(n=3, k=5, p=10, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @concept stiefel
#' @export
stiefel.runif <- function(n, k, p, type=c("list","array","riemdata")){
  ## PREPROCESSING
  N = round(n)
  p = round(p) # k-frame in R^p
  k = round(k)
  if (k > p){
    stop("* stiefel.runif : 'k <= p' is a required condition.")
  }
  ## MAIN COMPUTATION
  #  return object type
  retype = ifelse(missing(type), "riemdata",
                  match.arg(tolower(type), c("list","array","riemdata")))
  tmpout = runif_stiefel(p,k,N) # C++ version

  ## RETURN
  if (all(retype=="array")){
    return(tmpout)
  } else {
    tmpobj = wrap.stiefel(tmpout)
    if (all(retype=="list")){
      return(tmpobj$data)
    } else {
      return(tmpobj)
    }
  }
}

#  (02) stiefel.utest ==========================================================
#' Test of Uniformity on Stiefel Manifold
#' 
#' Given the data on Stiefel manifold \eqn{St(k,p)}, it tests whether the 
#' data is distributed uniformly.
#' 
#' @param stobj a S3 \code{"riemdata"} class for \eqn{N} Stiefel-valued data.
#' @param method (case-insensitive) name of the test method containing \describe{
#' \item{\code{"Rayleigh"}}{original Rayleigh statistic.}
#' \item{\code{"RayleighM"}}{modified Rayleigh statistic with better order of error.}
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
#' @seealso \code{\link{wrap.stiefel}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #   Compare Rayleigh's original and modified versions of the test
#' # 
#' # Test 1. sample uniformly from St(2,4)
#' # Test 2. use perturbed principal components from 'iris' data in R^4
#' #         which is concentrated around a point to reject H0.
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' #  1. uniform data
#' myobj1 = stiefel.runif(n=100, k=2, p=4)
#' 
#' #  2. perturbed principal components
#' data(iris)
#' irdat = list()
#' for (n in 1:100){
#'    tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
#'    irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
#' }
#' myobj2 = wrap.stiefel(irdat)
#' 
#' ## TEST
#' #  1. uniform data
#' stiefel.utest(myobj1, method="Rayleigh")
#' stiefel.utest(myobj1, method="RayleighM")
#' 
#' #  2. concentrated data
#' stiefel.utest(myobj2, method="rayleIgh")   # method names are 
#' stiefel.utest(myobj2, method="raYleiGhM")  # CASE - INSENSITIVE !
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @concept stiefel
#' @export
stiefel.utest <- function(stobj, method=c("Rayleigh","RayleighM")){
  ## CHECK INPUT
  FNAME      = "stiefel.utest"
  DNAME      = deparse(substitute(stobj)) # borrowed from HDtest
  check_inputmfd(stobj, FNAME)
  method.now = ifelse(missing(method),"rayleigh",
                      match.arg(tolower(method),c("rayleigh","rayleighm")))
  
  ## COMPUTATION
  #  prepare data in 3d array
  data3d <- aux_rmat2array3d(stobj)
  output <- switch(method.now,
                   rayleigh  = st.utest.Rayleigh(data3d, DNAME, is.modified = FALSE),
                   rayleighm = st.utest.Rayleigh(data3d, DNAME, is.modified = TRUE))
  return(output)
}
#' @keywords internal
#' @noRd
st.utest.Rayleigh <- function(x, dname, is.modified=FALSE){
  # Take the version from RiemStiefel {will be deprecated}
  p = dim(x)[1]      # r-frame in R^p with n observations
  r = dim(x)[2]
  n = dim(x)[3]
  xbar = array(0,c(p,r))
  for (i in 1:p){
    for (j in 1:r){
      xbar[i,j] = base::mean(as.vector(x[i,j,]))
    }
  }
  S    = p*n*sum(diag(t(xbar)%*%xbar))
  
  if (is.modified){
    hname   = "Modified Rayleigh Test of Uniformity on Stiefel Manifold"
    term1   = 1/(2*n)
    term2   = 1 - (S/((p*r) + 2))
    thestat = S*(1 - term1*term2)
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)  
  } else {
    hname   = "Rayleigh Test of Uniformity on Stiefel Manifold"
    thestat = S
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  }
  
  # COMPUTATION : DETERMINATION
  Ha      = paste("data is not uniformly distributed on St(",r,",",p,").",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}


#  (03) stiefel.optSA ==========================================================
#' Simulated Annealing on Stiefel Manifold
#' 
#' Simulated Annealing is a black-box, derivative-free optimization algorithm 
#' that iterates via stochastic search in the neighborhood of current position. 
#' \code{stiefel.optSA} solves the following problem
#' \deqn{\min_X f(X),\quad X \in St(p,k)}
#' without any other auxiliary information such as gradient or hessian involved. 
#' 
#' @param func a function to be \emph{minimized}.
#' @param p dimension parameter as in \eqn{St(k,p)}.
#' @param k dimension parameter as in \eqn{St(k,p)}.
#' @param ... extra parameters for SA algorithm including\describe{
#' \item{n.start}{number of runs; algorithm is executed \code{n.start} times (default: 5).}
#' \item{stepsize}{size of random walk on each component (default: 0.1).}
#' \item{maxiter}{maximum number of iterations for each run (default: 100).}
#' \item{cooling}{triplet for cooling schedule. See the section for the usage.}
#' \item{init.val}{if \code{NULL}, starts from a random point. Otherwise, a Stiefel matrix of size \eqn{(p,k)} should be provided for fixed starting point.}
#' \item{print.progress}{a logical; if \code{TRUE}, it prints each iteration.}
#' }
#' 
#' 
#' @return a named list containing: \describe{
#' \item{cost}{minimized function value.}
#' \item{solution}{a \eqn{(p\times k)} matrix that attains the \code{cost}.}
#' \item{accfreq}{frequency of acceptance moves.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #               Optimization for Eigen-Decomposition
#' #
#' # Given (5x5) covariance matrix S, eigendecomposition is indeed 
#' # an optimization problem cast on the stiefel manifold. Here, 
#' # we are trying to find top 3 eigenvalues and compare.
#' #-------------------------------------------------------------------
#' ## PREPARE
#' set.seed(121)                         # set seed
#' A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
#' myfunc <- function(p){                # cost function to minimize
#'   return(sum(-diag(t(p)%*%A%*%p)))
#' } 
#' 
#' ## SOLVE THE OPTIMIZATION PROBLEM
#' Aout = stiefel.optSA(myfunc, p=5, k=3, n.start=40, maxiter=200)
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
#' legend(1, 1, legend=c("optimization","decomposition"), col=c("red","blue"),
#'        lty=rep(1,2), pch=19)
#' par(opar)
#' 
#' @concept stiefel
#' @export
stiefel.optSA <- function(func, p, k, ...){
  # Preprocessing
  # 1. function
  if (!is.function(func)){
    stop("* stiefel.optSA : an input 'func' should be a function.")
  }
  # 2. implicit parameters
  params = list(...)
  pnames = names(params)
  
  n.start  = ifelse(("n.start"%in%pnames), max(1, params$n.start), 5) 
  stepsize = ifelse(("stepsize"%in%pnames), max(sqrt(.Machine$double.eps), 
                                                  params$stepsize), 0.1)
  maxiter  = ifelse(("maxiter"%in%pnames), max(10, round(params$maxiter)), 100)
  print.progress = ifelse(("print.progress"%in%pnames), as.logical(params$print.progress), FALSE)
  
  if ("init.val" %in% pnames){
    init.val = params$init.val
  } else {
    init.val = NULL
  }
  if ("cooling" %in% pnames){
    cooling = params$cooling
  } else {
    cooling =  c("exponential",10,0.9)
  }

  # 3. size
  size = c(round(p), round(k))
  if ((length(init.val)==0)&&(is.null(init.val))){ # not given
    if ((!is.vector(size))||(length(size)!=2)){
      stop("* stiefel.optSA : 'size' should be a vector of length 2.") 
    }  
    # for Stiefel Manifold Only
    n = round(size[1])
    p = round(size[2])
    initflag = FALSE
  } else { # if given, override size information
    n = nrow(init.val)
    p = ncol(init.val)
    initflag = TRUE 
  }
  # 3. other parameters
  my.nstart      = round(n.start)
  my.stepsize    = as.double(stepsize)
  my.temperature = sa_check_cooling(cooling, maxiter)
  
  ##############################################################
  # Main Run
  if (isTRUE(initflag)){
    out.now = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
    if (print.progress){
      print(paste0("* stiefel.optSA : iteration 1/",n.start," complete.."))
    }
  } else {
    init.val = base::qr.Q(base::qr(matrix(stats::rnorm(n*p),nrow=n)))
    out.now  = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
    if (print.progress){
      print(paste0("* stiefel.optSA : iteration 1/",n.start," complete.."))
    }
  }
  if (my.nstart > 1){
    for (it in 1:(my.nstart-1)){
      if (isTRUE(initflag)){
        out.tmp = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
      } else {
        init.val = base::qr.Q(base::qr(matrix(stats::rnorm(n*p),nrow=n)))
        out.tmp  = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
      }
      if (out.tmp$cost <= out.now$cost){ # update with a better one
        out.now = out.tmp
      }
      if (print.progress){
        print(paste0("* stiefel.optSA : iteration ",it+1,"/",n.start," complete.."))
      }
    }
  }
  
  ##############################################################
  # Report
  return(out.now)
}
# SA functions ------------------------------------------------------------
# (1) sa_check_cooling  : check the cooling schedule
#                         c("exponential", C, alpha)
#                         c("logarithmic", C, D)
#                         c("turnpike",    C, D)
# (2) sa_engine_Stiefel : working version of the function for Stiefel Manifold
#' @keywords internal
#' @noRd
sa_check_cooling <- function(coolschedule, itermax){
  if ((!is.vector(coolschedule))||(length(coolschedule)!=3)){
    stop("* stiefel.optSA : 'cooling' schedule must be a vector of length 3.")
  }
  themethod = coolschedule[1]
  itermax   = round(itermax)
  C         = as.double(coolschedule[2])
  if (C <= 0){
    stop("* stiefel.optSA : 'C' should be a nonnegative number.")
  }
  if (all(tolower(themethod)=="exponential")){
    alpha = as.double(coolschedule[3])
    if ((alpha <=0)||(alpha >=1)){
      stop("* stiefel.optSA : when the cooling schedule is 'exponential', 'alpha' should be a value in (0,1).")
    }
    outvec = C*(alpha^(1:itermax))
  } else if (all(tolower(themethod)=="logarithmic")){
    D = as.double(coolschedule[3])
    if (D <= 0){
      stop("* stiefel.optSA : when the cooling schedule is 'logarithmic', 'D' should be a positive real number.")
    }
    outvec = C/(log((1:itermax)+D))
  } else if (all(tolower(themethod)=="turnpike")){
    D = as.double(coolschedule[3])
    if (D <= 1){
      stop(" stiefel.optSA : when the cooling schedule is 'turnpike', 'D' should be a real number larger than 1.")
    }
    outvec = C*((D-1)/(log((1:itermax))))
  }
  
  if (any(is.infinite(outvec))){
    outvec[is.infinite(outvec)] = 2*max(outvec[!is.infinite(outvec)])
  }
  return(outvec)
}

# (2) sa_engine_Stiefel : working version of the function for Stiefel Manifold
#' @keywords internal
#' @noRd
sa_engine_Stiefel <- function(func, init.mat, temparature, stepsize){
  # initialization
  n = nrow(init.mat)
  p = ncol(init.mat)
  Eold    = func(init.mat)
  maxiter = length(temparature)
  
  # iteration
  sol.old = init.mat
  count   = 0
  for (k in 1:maxiter){
    # 1. generation unique to Stiefel 
    sol.tmp = sol.old + matrix(stats::rnorm(n*p, sd=stepsize), nrow=n)
    sol.tmp = base::qr.Q(base::qr(sol.tmp))
    # 2. energy evaluation & current temperature
    Etmp = func(sol.tmp)
    # 3. decision branching
    if (Etmp <= Eold){ # unconditional accept
      Enew    = Etmp
      sol.new = sol.tmp
      count   = count + 1
    } else {           # conditional accept
      Tk = temparature[k]
      tprob = exp((-(Etmp-Eold))/Tk)
      if (as.double(stats::runif(1)) < tprob){
        Enew    = Etmp
        sol.new = sol.tmp
        count   = count + 1
      } else {
        Enew    = Eold
        sol.new = sol.old
      }
    }
    # 4. update
    Eold    = Enew
    sol.old = sol.new
  }
  
  # report the run
  output = list()
  output$cost     = Eold
  output$solution = sol.old
  output$accfreq  = (count/maxiter)
  return(output)
}
