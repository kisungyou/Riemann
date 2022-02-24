## SOME SPECIAL FUNCTIONS ON SPHERE
#  (01) sphere.runif
#  (02) sphere.utest
#  (03) sphere.convert
#
#
#  (91) mixspnorm    & predict method
#  (92) mixsplaplace & predict method

# (01) sphere.runif ============================================================
#' Generate Uniform Samples on Sphere
#' 
#' It generates \eqn{n} random samples from \eqn{\mathcal{S}^{p-1}}. For convenient 
#' usage of users, we provide a number of options in terms of the return type.
#' 
#' @param n number of samples to be generated.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of length-\eqn{p} vectors.}
#' \item{\code{"matrix"}}{a \eqn{(n\times p)} where rows are unit vectors.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.sphere}} for more details (\emph{Default}).}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                       Draw Samples on Sphere
#' #
#' # Multiple return types on S^4 in R^5
#' #-------------------------------------------------------------------
#' dat.list = sphere.runif(n=10, p=5, type="list")
#' dat.matx = sphere.runif(n=10, p=5, type="matrix")
#' dat.riem = sphere.runif(n=10, p=5, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept sphere
#' @export
sphere.runif <- function(n, p, type=c("list","matrix","riemdata")){
  # PREPROCESSING
  # parameters
  N = round(n)
  p = round(p) # S^{p-1} in R^p
  
  # return object type
  if (missing(type)){
    retype = "riemdata"
  } else {
    retype = match.arg(tolower(type), c("list","matrix","riemdata"))
  }
  
  # GENERATION
  outmat = runif_sphere(N,p) 
  if (all(retype=="matrix")){
    return(outmat)
  } else {
    outobj = wrap.sphere(outmat)
    if (all(retype=="list")){
      return(outobj$data)
    } else {
      return(outobj)
    }
  }
}

# (02) sphere.utest ============================================================
#' Test of Uniformity on Sphere
#' 
#' Given \eqn{N} observations \eqn{\lbrace X_1, X_2, \ldots, X_M \brace} on 
#' \eqn{\mathcal{S}^{p-1}}, it tests whether the data is distributed uniformly 
#' on the sphere. 
#' 
#' @param spobj a S3 \code{"riemdata"} class for \eqn{N} Sphere-valued data.
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
#' 
#' @examples
#' #-------------------------------------------------------------------
#' #   Compare Rayleigh's original and modified versions of the test
#' #-------------------------------------------------------------------
#' #  Data Generation
#' myobj = sphere.runif(n=100, p=5, type="riemdata")
#' 
#' #  Compare 2 versions : Original vs Modified Rayleigh
#' sphere.utest(myobj, method="rayleigh")
#' sphere.utest(myobj, method="rayleighm")
#' 
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept sphere
#' @export
sphere.utest <- function(spobj, method=c("Rayleigh","RayleighM")){
  ## CHECK INPUT
  check_inputmfd(spobj, "sphere.utest")
  DNAME      = deparse(substitute(spobj)) # borrowed from HDtest
  method.now = ifelse(missing(method),"rayleigh",
                      match.arg(tolower(method),
                                c("rayleigh","rayleighm")))
  
  ## COMPUTE
  # prepare data in 3d array
  data3d <- aux_rvec2array3d(spobj)
  output <- switch(method.now,
                   rayleigh  = sp.utest.Rayleigh(data3d, DNAME, is.original = TRUE),
                   rayleighm = sp.utest.Rayleigh(data3d, DNAME, is.original = FALSE))
  return(output)
}
#' @keywords internal
#' @noRd
sp.utest.Rayleigh <- function(x, dname, is.original=TRUE){
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
  
  if (is.original){
    hname   = "Rayleigh Test of Uniformity on Sphere"
    thestat = S
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  } else {
    hname   = "Modified Rayleigh Test of Uniformity on Sphere"
    term1   = 1/(2*n)
    term2   = 1 - (S/((p*r) + 2))
    thestat = S*(1 - term1*term2)
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  }
  
  # COMPUTATION : DETERMINATION
  Ha      = paste("data is not uniformly distributed on ",p-1,"-sphere.",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}

#  (03) sphere.convert =========================================================
#  https://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates
#' Convert between Cartesian Coordinates and Geographic Coordinates
#' 
#' In geospatial data analysis, it is common to consider locations on the Earth as 
#' data. These locations, usually provided by latitude and longitude, are not directly 
#' applicable for spherical data analysis. We provide two functions - \code{sphere.geo2xyz} and \code{sphere.xyz2geo} - 
#' that convert geographic coordinates in longitude/latitude into a unit-norm vector on \eqn{\mathcal{S}^2}, and vice versa. 
#' As a convention, latitude and longitude are represented as \emph{decimal degrees}. 
#' 
#' @param lat latitude (in decimal degrees).
#' @param lon longitude (in decimal degrees).
#' @param xyz a unit-norm vector in \eqn{\mathcal{S}^{2}}.
#' 
#' @return transformed data.
#' 
#' @examples 
#' ## EXAMPLE DATA WITH POPULATED US CITIES
#' data(cities)
#' 
#' ## SELECT ALBUQUERQUE
#' geo = cities$coord[1,]
#' xyz = cities$cartesian[1,]
#' 
#' ## CHECK TWO INPUT TYPES AND THEIR CONVERSIONS
#' sphere.geo2xyz(geo[1], geo[2])
#' sphere.xyz2geo(xyz)
#' 
#' @name sphere.convert
#' @concept sphere
#' @rdname sphere.convert
NULL

#' @rdname sphere.convert
#' @export
sphere.geo2xyz <- function(lat, lon){
  xlat = as.double(lat)*pi/180
  xlon = as.double(lon)*pi/180
  
  x = cos(xlat)*cos(xlon)
  y = cos(xlat)*sin(xlon)
  z = sin(xlat)
  
  outvec = c(x,y,z)
  return(outvec/sqrt(sum(outvec^2)))
}

#' @rdname sphere.convert
#' @export
sphere.xyz2geo <- function(xyz){
  x = as.double(xyz[1])
  y = as.double(xyz[2])
  z = as.double(xyz[3])
  
  latitude  = 180*asin(z)/pi
  longitude = 180*atan2(y,x)/pi
  
  output = c(latitude, longitude)
  return(output)
}





# (92) mixsplaplace & predict method --------------------------------------
#' Finite Mixture of Spherical Laplace Distributions
#' 
#' For \eqn{n} observations on a \eqn{(p-1)} sphere in \eqn{\mathbf{R}^p}, 
#' a finite mixture model is fitted whose components are spherical Laplace distributions via the following model
#' \deqn{f(x; \left\lbrace w_k, \mu_k, \sigma_k \right\rbrace_{k=1}^K) = \sum_{k=1}^K w_k SL(x; \mu_k, \sigma_k)}
#' with parameters \eqn{w_k}'s for component weights, \eqn{\mu_k}'s for component locations, and \eqn{\sigma_k}'s for component scales. 
#' 
#' @param data data vectors in form of either an \eqn{(n\times p)} matrix or a length-\eqn{n} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param k the number of clusters (default: 2).
#' @param same.sigma a logical; \code{TRUE} to use same scale parameter across all components, or \code{FALSE} otherwise.
#' @param variants type of the class assignment methods, one of \code{"soft"},\code{"hard"}, and \code{"stochastic"}.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 50).}
#' \item{eps}{stopping criterion for the EM algorithm (default: 1e-6).}
#' \item{printer}{a logical; \code{TRUE} to show history of the algorithm, \code{FALSE} otherwise.}
#' }
#' 
#' @return a named list of S3 class \code{mixsplaplace} containing
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
#' \item{loglkd}{log likelihood of the fitted model.}
#' \item{criteria}{a vector of information criteria.}
#' \item{parameters}{a list containing \code{proportion}, \code{location}, and \code{scale}. See the section for more details.}
#' \item{membership}{an \eqn{(n\times k)} row-stochastic matrix of membership.}
#' }
#' 
#' @examples 
#' \donttest{
#' # LOAD THE CITY DATA AND WRAP AS RIEMOBJ
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # FIT THE MODEL WITH DIFFERENT K's
#' k2 = mixsplaplace(locations, k=2)
#' k3 = mixsplaplace(locations, k=3)
#' k4 = mixsplaplace(locations, k=4)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed2, col=k2$cluster, pch=19, main="K=2")
#' plot(embed2, col=k3$cluster, pch=19, main="K=3")
#' plot(embed2, col=k4$cluster, pch=19, main="K=4")
#' par(opar)
#' }
#' 
#' @section Parameters of the fitted model:
#' A fitted model is characterized by three parameters. For \eqn{k}-mixture model on a \eqn{(p-1)} 
#' sphere in \eqn{\mathbf{R}^p}, (1) \code{proportion} is a length-\eqn{k} vector of component weight 
#' that sums to 1, (2) \code{location} is an \eqn{(k\times p)} matrix whose rows are per-cluster locations, and 
#' (3) \code{concentration} is a length-\eqn{k} vector of scale parameters for each component.
#' 
#' @concept sphere
#' @export
mixsplaplace <- function(data, k=2, same.sigma=FALSE, variants=c("soft","hard","stochastic"), ...){
  ## PREPROCESSING -------------------------------------------------------------
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  FNAME  = "mixspnorm"
  
  pars   = list(...)
  pnames = names(pars)
  
  if ("maxiter"%in% pnames){
    myiter = max(pars$maxiter, 50)
  } else {
    myiter = 100
  }
  if ("eps" %in% pnames){
    myeps = max(.Machine$double.eps, as.double(pars$eps))
  } else {
    myeps = 1e-6
  }
  if ("printer" %in% pnames){
    myprint = as.logical(pars$printer)
  } else {
    myprint = FALSE
  }
  same.sigma = as.logical(same.sigma)
  
  myn = base::nrow(x)
  myp = base::ncol(x)-1
  myk = max(1, round(k))
  myvers = match.arg(variants)
  
  ## INITIALIZATION ------------------------------------------------------------
  #  label
  initlabel = as.vector(stats::kmeans(x, myk, nstart=round(5))$cluster)
  #  membership
  old.eta = array(0,c(myn,myk))
  for (i in 1:myn){
    old.eta[i,initlabel[i]] = 1 # {0,1} for the initial
  }
  old.mu     = moSL.median(spobj, old.eta)
  old.mobj   = wrap.sphere(old.mu)
  old.d2med  = moSL.d2medmat(spobj, old.mobj)
  old.sigma  = moSL.updatesig(old.eta, old.d2med, myp, homogeneous = same.sigma)
  old.pi     = as.vector(base::colSums(old.eta)/myn)
  old.loglkd = moSL.loglkd(old.d2med, old.sigma, old.pi, myp)
  
  ## ITERATION -----------------------------------------------------------------
  inc.params = rep(0, 5)
  for (it in 1:myiter){
    # E-step
    new.eta = moSL.eta(old.d2med, old.sigma, old.pi, myp)
    
    # H/S-step by option
    if (all(myvers=="hard")){
      new.eta = spmix.hard(new.eta)
    } else if (all(myvers=="stochastic")){
      new.eta = spmix.stochastic(new.eta)
    }
    
    # Stop if there is empty cluster
    new.label = apply(new.eta, 1, which.max)
    if (length(unique(new.label)) < myk){
      break
    }
    
    # M-step
    # M1. mu / centers & d2med
    new.mu    = moSL.median(spobj, new.eta)
    new.mobj  = wrap.sphere(new.mu)
    new.d2med = moSL.d2medmat(spobj, new.mobj)
    # M2. sigma / scales
    new.sigma = moSL.updatesig(new.eta, new.d2med, myp, homogeneous = same.sigma)
    # M3. proportions
    new.pi    = as.vector(base::colSums(new.eta)/myn)
    # update
    new.loglkd = moSL.loglkd(new.d2med, new.sigma, new.pi, myp)
    
    # Incremental changes
    inc.params[1] = base::norm(old.mu-new.mu, type = "F")
    inc.params[2] = base::sqrt(base::sum((old.sigma - new.sigma)^2))
    inc.params[3] = base::sqrt(base::sum((old.pi-new.pi)^2))
    inc.params[4] = base::norm(old.d2med-new.d2med, type="F")
    inc.params[5] = base::norm(old.eta-new.eta, type="F")
    
    # rule : log-likelihood should increase
    if (new.loglkd < old.loglkd){
      if (myprint){
        print(paste0("* mixsplaplace : terminated at iteration ", it, " : log-likelihood is decreasing."))
      }
      break
    } else {
      old.eta    = new.eta
      old.mu     = new.mu
      old.mobj   = new.mobj
      old.d2med  = new.d2med
      old.sigma  = new.sigma
      old.pi     = new.pi
      old.loglkd = new.loglkd
    }
    if (max(inc.params) < myeps){
      if (myprint){
        print(paste0("* mixsplaplace : terminated at iteration ", it," : all parameters converged."))  
      }
      break
    }
    if (myprint){
      print(paste0("* mixsplaplace : iteration ",it,"/",myiter," complete with loglkd=",round(old.loglkd,5),"."))  
    }
  }
  
  ## INFORMATION CRITERION -----------------------------------------------------
  if (!same.sigma){
    par.k = myp*myk + myk + (myk-1)
  } else {
    par.k = myp*myk + 1 + (myk-1)
  }
  
  AIC = -2*old.loglkd + 2*par.k
  BIC = -2*old.loglkd + par.k*log(myn)
  HQIC = -2*old.loglkd + 2*par.k*log(log(myn))
  AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
  
  infov = c(AIC, AICc, BIC, HQIC)
  names(infov) = c("AIC","AICc","BIC","HQIC")
  
  ## RETURN --------------------------------------------------------------------
  output = list()
  output$cluster  = spmix.getcluster(old.eta)
  output$loglkd   = old.loglkd
  output$criteria = infov
  output$parameters = list(proportion=old.pi, location=old.mu, scale=old.sigma)
  output$membership = old.eta
  return(structure(output, class="mixsplaplace"))
}



# auxiliary functions -----------------------------------------------------
# 1. moSL.median    : weighted Frechet medians
# 2. moSL.d2medmat  : compute data-median pairwise distance
# 3. moSL.solvesig  : minimize C/sigma + log(C(sigma))
# 4. moSL.updatesig : update sigma (scale) parameters
# 5. moSL.loglkd    : log-likelihood
# 6. moSL.eta       : compute soft membership

# 1. compute weighted Frechet median : deals with Riemdata object
#' @keywords internal
#' @noRd
moSL.median <- function(spobj, membership){
  N = base::length(spobj$data)
  P = base::length(as.vector(spobj$data[[1]]))
  K = base::ncol(membership)
  
  output = array(0,c(K,P))
  for (k in 1:K){
    partweight = as.vector(membership[,k])
    if (sum(partweight==1)==1){
      output[k,] = as.vector(spobj$data[[which(partweight==1)]])
    } else {
      sel_id     = (partweight > .Machine$double.eps)
      sel_data   = spobj$data[sel_id]
      sel_weight = partweight[sel_id]
      sel_weight = sel_weight/base::sum(sel_weight)
      output[k,] = as.vector(inference_median_intrinsic("sphere", sel_data, sel_weight, 100, 1e-6)$median)
    }
  }
  return(output)
}
# 2. moSL.d2medmat 
#' @keywords internal
#' @noRd
moSL.d2medmat <- function(obj.data, obj.medians){
  d2sqmat = as.matrix(basic_pdist2("sphere", obj.data$data, obj.medians$data, "intrinsic"))
  return(d2sqmat) 
}
# 3. moSL.solvesig : minimize C/sigma + log(C(sigma))
#' @keywords internal
#' @noRd
moSL.solvesig <- function(C, p){
  # cost function
  fun_cost <- function(sigma){
    # term : first
    out1 = C/sigma
    
    # term : second
    myfunc <- function(r){
      return(exp(-r/sigma)*(sin(r)^(p-1)))
    }
    t1   = 2*(pi^(p/2))/(gamma(p/2))
    t2   = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    out2 = log(t1) + log(t2)
    
    # return the output
    return(out1+out2)
  }
  
  # minimization
  myint  = c(0.01, 100)*C
  output = as.double(stats::optimize(fun_cost, interval=myint, maximum=FALSE, tol=1e-6)$minimum)
  return(output)
}
# 4. moSL.updatesig : update sigma (scale) parameters
#' @keywords internal
#' @noRd
moSL.updatesig <- function(membership, d2mat, p, homogeneous=TRUE){
  N = base::nrow(membership)
  K = base::ncol(membership)
  mysigs = rep(0, K)
  
  if (homogeneous){ # homogeneous
    A = base::sum(d2mat*membership)
    B = base::sum(membership)
    C = (A/B)
    sig.single = moSL.solvesig(C, p)
    for (k in 1:K){
      mysigs[k] = sig.single
    }
  } else {          # heterogeneous
    for (k in 1:K){
      A = sum(as.vector(d2mat[,k])*as.vector(membership[,k]))
      B = sum(as.vector(membership[,k]))
      C = (A/B)
      
      mysigs[k] = moSL.solvesig(C, p)
    }
  }
  return(mysigs)
}
# 5. moSL.loglkd    : log-likelihood
#' @keywords internal
#' @noRd
moSL.loglkd <- function(d2med, sigmas, props, p){
  N = base::nrow(d2med)
  K = base::ncol(d2med)
  
  # function to evaluate normalizing constant
  eval_constant <- function(sigma){
    myfunc <- function(r){
      return(exp(-r/sigma)*(sin(r)^(p-1)))
    }
    t1   = 2*(pi^(p/2))/(gamma(p/2))
    t2   = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    return(t1*t2) 
  }
  
  # normalizing constants per class
  vecCsig = rep(0, K)
  for (k in 1:K){
    vecCsig[k] = eval_constant(sigmas[k])
  }
  
  # evaluate the density
  mixdensity = array(0,c(N,K))
  for (n in 1:N){
    for (k in 1:K){
      mixdensity[n,k] = props[k]*exp(-(d2med[n,k])/sigmas[k])/vecCsig[k]
    }
  }
  
  # evaluate the output
  return(base::sum(base::log(base::rowSums(mixdensity))))
}
# 6. moSL.eta 
#' @keywords internal
#' @noRd
moSL.eta <- function(d2med, sigmas, props, p){
  N = base::nrow(d2med)
  K = base::ncol(d2med)
  
  # function to evaluate normalizing constant
  eval_constant <- function(sigma){
    myfunc <- function(r){
      return(exp(-r/sigma)*(sin(r)^(p-1)))
    }
    t1   = 2*(pi^(p/2))/(gamma(p/2))
    t2   = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    return(t1*t2) 
  }
  
  # normalizing constants per class
  vecCsig = rep(0, K)
  for (k in 1:K){
    vecCsig[k] = eval_constant(sigmas[k])
  }
  
  # evaluate
  output = array(0,c(N,K))
  for (k in 1:K){
    tgtdvec    = as.vector(d2med[,k])
    output[,k] = exp((-tgtdvec/sigmas[k]) - base::log(vecCsig[k]) + base::log(props[k]))
  }
  for (n in 1:N){
    tgtrow = as.vector(output[n,])
    output[n,] = tgtrow/base::sum(tgtrow)
  }
  return(output)
}



#' S3 Method for Prediction upon Fitted 'mixsplaplace' Model
#' 
#' Given new data with the fitted mixture of spherical Laplace distributions on a sphere, predict the 
#' class labels for the newly provided data according to the fitted model. 
#' 
#' @param object an object of \code{mixsplaplace} class. See \code{\link{mixsplaplace}} for more details.
#' @param newdata data vectors in form of either an \eqn{(m\times p)} matrix or a length-\eqn{m} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return a length-\eqn{m} vector of class labels.
#' 
#' @examples 
#' \donttest{
#' # LOAD THE CITY DATA AND WRAP AS RIEMOBJ
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#'  }
#'  
#' # FIT THE MODEL K=3
#' k3fit     = mixsplaplace(locations, k=3)
#' k3fitlab  = k3fit$cluster
#' 
#' # PREDICT THE CLASS LABEL WITH THE SAME DATA
#' k3predict = predict(k3fit, locations)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(embed2, col=k3fitlab,  pch=19, main="fitted")
#' plot(embed2, col=k3predict, pch=19, main="predicted")
#' par(opar)
#' }
#' 
#' @seealso \code{\link{mixsplaplace}}
#' @concept sphere
#' @export
predict.mixsplaplace <- function(object, newdata, ...){
  # PREPARE
  if (!inherits(object, "mixsplaplace")){
    stop("* predict : input is not an object of 'mixsplaplace' class.")
  }
  spobj = wrap.sphere(newdata)
  x     = sp2mat(spobj)
  myp   = base::ncol(x)-1
  
  # COMPUTE : Cluster Label Only
  old.mu    = object$parameters$location
  old.mobj  = wrap.sphere(old.mu)
  old.sigma = object$parameters$scale
  old.pi    = object$parameters$proportion
  old.d2med = moSL.d2medmat(spobj, old.mobj)
  
  fin.eta   = moSL.eta(old.d2med, old.sigma, old.pi, myp)
  output    = spmix.getcluster(fin.eta)
  return(output)
}