## SOME SPECIAL FUNCTIONS ON SPHERE
#  (01) sphere.runif
#  (02) sphere.utest
#  (03) sphere.convert
#
#
#  (91) mixspnorm & predict method

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






# (91) mixspnorm & predict method -----------------------------------------
#' Finite Mixture of Spherical Normal Distributions
#' 
#' For \eqn{n} observations on a \eqn{(p-1)} sphere in \eqn{\mathbf{R}^p}, 
#' a finite mixture model is fitted whose components are spherical normal distributions via the following model
#' \deqn{f(x; \left\lbrace w_k, \mu_k, \lambda_k \right\rbrace_{k=1}^K) = \sum_{k=1}^K w_k SN(x; \mu_k, \lambda_k)}
#' with parameters \eqn{w_k}'s for component weights, \eqn{\mu_k}'s for component locations, and \eqn{\lambda_k}'s for component concentrations. 
#' 
#' @param data data vectors in form of either an \eqn{(n\times p)} matrix or a length-\eqn{n} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param k the number of clusters (default: 2).
#' @param same.lambda a logical; \code{TRUE} to use same concentration parameter across all components, or \code{FALSE} otherwise.
#' @param variants type of the class assignment methods, one of \code{"soft"},\code{"hard"}, and \code{"stochastic"}.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 50).}
#' \item{eps}{stopping criterion for the EM algorithm (default: 1e-6).}
#' \item{printer}{a logical; \code{TRUE} to show history of the algorithm, \code{FALSE} otherwise.}
#' }
#' 
#' @return a named list of S3 class \code{mixspnorm} containing
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
#' \item{loglkd}{log likelihood of the fitted model.}
#' \item{criteria}{a vector of information criteria.}
#' \item{parameters}{a list containing \code{proportion}, \code{center}, and \code{concentration}. See the section for more details.}
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
#' k2 = mixspnorm(locations, k=2)
#' k3 = mixspnorm(locations, k=3)
#' k4 = mixspnorm(locations, k=4)
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
#' that sums to 1, (2) \code{center} is an \eqn{(k\times p)} matrix whose rows are cluster centers, and 
#' (3) \code{concentration} is a length-\eqn{k} vector of concentration parameters for each component.
#' 
#' @concept sphere
#' @export
mixspnorm <- function(data, k=2, same.lambda=FALSE, variants=c("soft","hard","stochastic"), ...){
  ## PREPROCESSING
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  FNAME  = "mixspnorm"
  
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 100))
  myeps  = min(1e-6, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-6)))
  myprint = ifelse(("printer"%in%pnames), as.logical(pars$printer), FALSE)
  same.lambda = as.logical(same.lambda)
  
  myn = nrow(x)
  myp = ncol(x)
  myk = max(round(k), 1)
  myvers = match.arg(variants)
  
  ## INITIALIZATION
  #  label
  initlabel = as.vector(stats::kmeans(x, myk, nstart=round(5))$cluster)
  #  membership
  old.eta = array(0,c(myn,myk))
  for (i in 1:myn){
    old.eta[i,initlabel[i]] = 1 # {0,1} for the initial
  }
  old.mu     = spmix.frechet(spobj, old.eta)
  old.mobj   = wrap.sphere(old.mu)
  old.d2mat  = spmix.d2sqmat(spobj, old.mobj)
  old.lambda = spmix.updatelbd(old.eta, old.d2mat, myp, homogeneous = same.lambda)
  old.pi     = as.vector(base::colSums(old.eta)/myn)
  old.loglkd = spmix.loglkd(old.d2mat, old.lambda, old.pi, myp)
  
  ## ITERATION
  inc.params = rep(0,5)
  for (it in 1:myiter){
    # E-Step
    new.eta = spmix.eta(old.d2mat, old.lambda, old.pi, myp)
    
    # # H/S-Step by Option
    if (all(myvers=="hard")){
      new.eta = spmix.hard(new.eta)
    } else if (all(myvers=="stochastic")){
      new.eta = spmix.stochastic(new.eta)
    }
    
    # M-Step
    # M1. mu / centers & d2mat
    new.mu     = spmix.frechet(spobj, new.eta) # (k,p) matrix of row-stacked centroids
    new.mobj   = wrap.sphere(new.mu)
    new.d2mat  = spmix.d2sqmat(spobj, new.mobj)
    # M2. lambda / concentration
    new.lambda = spmix.updatelbd(new.eta, new.d2mat, myp, homogeneous = same.lambda)
    # M3. proportion
    new.pi     = as.vector(base::colSums(new.eta)/myn)
    # update 
    new.loglkd = spmix.loglkd(new.d2mat, new.lambda, new.pi, myp)
    
    inc.params[1] = base::norm(old.mu-new.mu, type = "F")
    inc.params[2] = base::sqrt(base::sum((old.lambda - new.lambda)^2))
    inc.params[3] = base::sqrt(base::sum((old.pi-new.pi)^2))
    inc.params[4] = base::norm(old.d2mat-new.d2mat, type="F")
    inc.params[5] = base::norm(old.eta-new.eta, type="F")
    
    if (new.loglkd < old.loglkd){ # rule : loglkd should increase
      if (myprint){
        print(paste0("* mixspnorm : terminated at iteration ", it," : log-likelihood is decreasing."))  
      }
      break
    } else {
      old.eta    = new.eta
      old.mu     = new.mu
      old.mobj   = new.mobj
      old.d2mat  = new.d2mat
      old.lambda = new.lambda
      old.pi     = new.pi
      old.loglkd = new.loglkd
    }
    if (max(inc.params) < myeps){ # rule : all parameters converge -> stop
      if (myprint){
        print(paste0("* mixspnorm : terminated at iteration ", it," : all parameters converged."))  
      }
      break
    }
    if (myprint){
      print(paste0("* mixspnorm : iteration ",it,"/",myiter," complete with loglkd=",round(old.loglkd,5)))  
    }
  }
  
  ## INFORMATION CRITERION
  if (!same.lambda){
    par.k = ((myp-1)*myk) + myk + (myk-1)  
  } else {
    par.k = ((myp-1)*myk) + 1   + (myk-1)  
  }
  AIC = -2*old.loglkd + 2*par.k
  BIC = -2*old.loglkd + par.k*log(myn)
  HQIC = -2*old.loglkd + 2*par.k*log(log(myn))
  AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
  
  infov = c(AIC, AICc, BIC, HQIC)
  names(infov) = c("AIC","AICc","BIC","HQIC")
  
  ## RETURN
  output = list()
  output$cluster  = spmix.getcluster(old.eta)
  output$loglkd   = old.loglkd  # max loglkd
  output$criteria = infov       # min AIC/AICc/BIC/HQIC
  output$parameters = list(proportion=old.pi, center=old.mu, concentration=old.lambda)
  output$membership = old.eta
  return(structure(output, class="mixspnorm"))
}



# auxiliary renewed -------------------------------------------------------
#  1. spmix.frechet    : weighted frechet mean
#  2. spmix.d2sqmat    : compute data and mean pairwise distance squared
#  3. spmix.solvelbd   : solve min A*lbd + B*log(Z(lbd))
#  4. spmix.updatelbd  : update lambda parameters
#  5. spmix.loglkd     : compute log-likelihood
#  6. spmix.eta        : compute eta (membership)
#  7. spmix.getcluster : find the label of the data given the membership
#  8. spmix.hard & spmix.stochastic

#  1. compute weighted frechet mean : deals with Riemdata object
#' @keywords internal
#' @noRd
spmix.frechet <- function(spobj, membership){
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
      output[k,] = as.vector(inference_mean_intrinsic("sphere", sel_data, sel_weight, 100, 1e-6)$mean)
    }
  }
  return(output)
}
#  2. spmix.d2sqmat   : compute data and mean pairwise distance squared
#' @keywords internal
#' @noRd
spmix.d2sqmat <- function(obj.data, obj.mean){
  d2sqmat = as.matrix(basic_pdist2("sphere", obj.data$data, obj.mean$data, "intrinsic"))^2
  return(d2sqmat)
}
#  3. spmix.solvelbd : solve min A*lbd + B*log(Z(lbd)) in R^P
#' @keywords internal
#' @noRd
spmix.solvelbd <- function(A, B, P){ 
  myfun <- function(par.lbd){
    myintegral <- function(r){
      return(exp(-par.lbd*(r^2)/2)*((sin(r))^(P-2)))
    }
    t1 = 2*(pi^((P-1)/2))/gamma((P-1)/2) # one possible source of error
    t2 = stats::integrate(myintegral, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    logZlbd = base::log(t1*t2)
    output  = A*par.lbd + B*logZlbd
    return(output)
  }
  
  sol = stats::optimize(myfun, interval=c(1e-15, 1e+5), maximum = FALSE)$minimum
  return(as.double(sol))
}
#  4. spmix.updatelbd : update lambda parameters
#' @keywords internal
#' @noRd
spmix.updatelbd <- function(membership, d2mat, P, homogeneous=TRUE){
  N = base::nrow(membership)
  K = base::ncol(membership)
  mylbds = rep(0,K)
  
  if (homogeneous){
    A = base::sum(d2mat*membership)/2
    B = base::sum(membership)
    lbd.single = spmix.solvelbd(A, B, P)
    for (k in 1:K){
      mylbds[k] = lbd.single
    }
  } else {
    for (k in 1:K){
      A = base::sum(as.vector(d2mat[,k])*as.vector(membership[,k]))/2
      B = base::sum(as.vector(membership[,k]))
      mylbds[k] = spmix.solvelbd(A, B, P)
    }
  }
  return(mylbds)
}
#  5. spmix.loglkd    : compute log-likelihood
#' @keywords internal
#' @noRd
spmix.loglkd <- function(d2mat, lambdas, props, P){
  N = base::nrow(d2mat)
  K = base::ncol(d2mat)
  
  evalZ <- function(par.lbd){
    myintegral <- function(r){
      return(exp(-par.lbd*(r^2)/2)*((sin(r))^(P-2)))
    }
    lbdt1 = 2*(pi^((P-1)/2))/gamma((P-1)/2) # one possible source of error
    lbdt2 = stats::integrate(myintegral, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    Zlbd  = lbdt1*lbdt2
    return(Zlbd)
  }
  
  logZlbd = rep(0,K)
  logprop = base::log(props)
  for (k in 1:K){
    logZlbd[k] = base::log(evalZ(lambdas[k]))
  }
  
  term1 = base::sum(logprop)*N
  term2 = -base::sum(logZlbd)*N
  term3 = base::sum(d2mat%*%(-base::diag(lambdas)/2))
  return(term1+term2+term3)
}
# 6. spmix.eta : compute eta (membership)
#' @keywords internal
#' @noRd
spmix.eta <- function(d2mat, lambdas, props, P){
  N = base::nrow(d2mat)
  K = base::ncol(d2mat)
  
  evalZ <- function(par.lbd){
    myintegral <- function(r){
      return(exp(-par.lbd*(r^2)/2)*((sin(r))^(P-2)))
    }
    lbdt1 = 2*(pi^((P-1)/2))/gamma((P-1)/2) # one possible source of error
    lbdt2 = stats::integrate(myintegral, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    Zlbd  = lbdt1*lbdt2
    return(Zlbd)
  }
  
  vecZlbd = rep(0,K)
  for (k in 1:K){
    vecZlbd[k] = evalZ(lambdas[k])
  }
  
  output = array(0,c(N,K))
  for (k in 1:K){
    tgtd2vec   = as.vector(d2mat[,k])
    output[,k] = exp((-lambdas[k]*tgtd2vec/2) - base::log(vecZlbd[k]) + base::log(props[k]))
  }
  for (n in 1:N){
    tgtrow = as.vector(output[n,])
    output[n,] = tgtrow/base::sum(tgtrow)
  }
  return(output)
}
#  7. spmix.getcluster : find the label of the data given the membership
#' @keywords internal
#' @noRd
spmix.getcluster <- function(membership){
  N = base::nrow(membership)
  K = base::ncol(membership)
  
  output = rep(0,N)
  for (n in 1:N){
    tgtvec = as.vector(membership[n,])
    output[n] = base::which.max(tgtvec)
  }
  return(output)
}
#  8. spmix.hard & spmix.stochastic
#' @keywords internal
#' @noRd
spmix.hard <- function(membership){
  N = base::nrow(membership)
  K = base::ncol(membership)
  
  output = array(0,c(N,K))
  for (n in 1:N){
    tgt = as.vector(membership[n,])
    output[n,which.max(tgt)] = 1
  }
  return(output)
}
#' @keywords internal
#' @noRd
spmix.stochastic <- function(membership){
  N = base::nrow(membership)
  K = base::ncol(membership)
  
  output = array(0,c(N,K))
  vec1ks = (1:K)
  for (n in 1:N){
    probn = as.vector(membership[n,])
    output[n,base::sample(vec1ks, 1, prob=probn)] = 1
  }
  return(output)
}



## LOAD THE DATA AND WRAP AS RIEMOBJ
# data(cities)
# locations = cities$cartesian
# embed2    = array(0,c(60,2))
# for (i in 1:60){
#   embed2[i,] = sphere.xyz2geo(locations[i,])
# }
# 
# ## CLUSTERING
# k2int = mixspnorm(locations, k=2, same.lambda = TRUE, printer=TRUE)
# k3int = sphere.mixnorm(locations, k=3, same.lambda = TRUE, printer=TRUE)
# k5int = sphere.mixnorm(locations, k=5, same.lambda = TRUE, printer=TRUE)
# k2ext = predict(movMF::movMF(locations, k=2), locations)
# k3ext = predict(movMF::movMF(locations, k=3), locations)
# k5ext = predict(movMF::movMF(locations, k=5), locations)
# 
# par(mfrow=c(2,3))
# plot(embed2, col=k2int$cluster, pch=19)
# plot(embed2, col=k3int$cluster, pch=19)
# plot(embed2, col=k5int$cluster, pch=19)
# plot(embed2, col=k2ext, pch=19)
# plot(embed2, col=k3ext, pch=19)
# plot(embed2, col=k5ext, pch=19)

#' S3 Method for Prediction upon Fitted 'mixspnorm' Model
#' 
#' Given new data with the fitted mixture of spherical normals on a sphere, predict the 
#' class labels for the newly provided data according to the fitted model. 
#' 
#' @param object an object of \code{mixspnorm} class. See \code{\link{mixspnorm}} for more details.
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
#' k3fit     = mixspnorm(locations, k=3)
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
#' @seealso \code{\link{mixspnorm}}
#' @concept sphere
#' @export
predict.mixspnorm <- function(object, newdata, ...){
  # PREPARE
  if (!inherits(object, "mixspnorm")){
    stop("* predict : input is not an object of 'mixspnorm' class.")
  }
  spobj = wrap.sphere(newdata)
  x     = sp2mat(spobj)
  myp   = base::ncol(x)
  
  # COMPUTE : Cluster Label Only
  old.mu     = object$parameters$center
  old.mobj   = wrap.sphere(old.mu)
  old.lambda = object$parameters$concentration
  old.pi     = object$parameters$proportion
  old.d2mat  = spmix.d2sqmat(spobj, old.mobj) # (N x K)
  
  fin.eta = spmix.eta(old.d2mat, old.lambda, old.pi, myp)
  output  = spmix.getcluster(fin.eta)
  return(output)
}