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
#' @return a named list of S3 class \code{riemmix} containing
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
#' # ---------------------------------------------------- #
#' #                 FITTING THE MODEL
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit the model with different numbers of clusters
#' k2 = moSL(locations, k=2)
#' k3 = moSL(locations, k=3)
#' k4 = moSL(locations, k=4)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed2, col=k2$cluster, pch=19, main="K=2")
#' plot(embed2, col=k3$cluster, pch=19, main="K=3")
#' plot(embed2, col=k4$cluster, pch=19, main="K=4")
#' par(opar)
#' 
#' # ---------------------------------------------------- #
#' #                   USE S3 METHODS
#' # ---------------------------------------------------- #
#' # Use the same 'locations' data as new data 
#' # (1) log-likelihood
#' newloglkd = round(loglkd(k3, locations), 5)
#' fitloglkd = round(k3$loglkd, 5)
#' print(paste0("Log-likelihood for K=3 fitted    : ", fitloglkd))
#' print(paste0("Log-likelihood for K=3 predicted : ", newloglkd))
#' 
#' # (2) label
#' newlabel = label(k3, locations)
#' 
#' # (3) density
#' newdensity = density(k3, locations)
#' }
#' 
#' @section Parameters of the fitted model:
#' A fitted model is characterized by three parameters. For \eqn{k}-mixture model on a \eqn{(p-1)} 
#' sphere in \eqn{\mathbf{R}^p}, (1) \code{proportion} is a length-\eqn{k} vector of component weight 
#' that sums to 1, (2) \code{location} is an \eqn{(k\times p)} matrix whose rows are per-cluster locations, and 
#' (3) \code{concentration} is a length-\eqn{k} vector of scale parameters for each component.
#' 
#' @section Note on S3 methods:
#' There are three S3 methods; \code{loglkd}, \code{label}, and \code{density}. Given a random sample of 
#' size \eqn{m} as \code{newdata}, (1) \code{loglkd} returns a scalar value of the computed log-likelihood, 
#' (2) \code{label} returns a length-\eqn{m} vector of cluster assignments, and (3) \code{density} 
#' evaluates densities of every observation according ot the model fit. 
#' 
#' @concept sphere
#' @export
moSL <- function(data, k=2, same.sigma=FALSE, variants=c("soft","hard","stochastic"), ...){
  ## PREPROCESSING -------------------------------------------------------------
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  FNAME  = "moSL"
  
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
  return(structure(output, class=c("moSL","riemmix")))
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




# Methods -----------------------------------------------------------------
# () loglkd  : compute the log-likelihood
# () label   : predict the labels 
# () density : evaluate the density

#  S3 METHOD : LOGLKD
#' @param object a fitted \code{moSL} model from the \code{\link{moSL}} function.
#' @param newdata data vectors in form of either an \eqn{(m\times p)} matrix or a length-\eqn{m} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @rdname moSL
#' @concept sphere
#' @export
loglkd.moSL <- function(object, newdata){
  # PREPARE
  if (!inherits(object, "moSL")){
    stop("* loglkd.moSL : input is not an object of 'moSL' class.")
  }
  spobj = wrap.sphere(newdata)
  myp   = base::ncol(object$parameters$location)-1
  
  # INTERMEDIATE VALUES
  old.mu    = object$parameters$location
  old.mobj  = wrap.sphere(old.mu)
  old.sigma = object$parameters$scale
  old.pi    = object$parameters$proportion
  old.d2med = moSL.d2medmat(spobj, old.mobj)
  
  # COMPUTE AND RETURN
  return(moSL.loglkd(old.d2med, old.sigma, old.pi, myp))
}

# S3 METHOD : LABEL
#' @rdname moSL
#' @concept sphere
#' @export
label.moSL <- function(object, newdata){
  # PREPARE
  if (!inherits(object, "moSL")){
    stop("* label.moSL : input is not an object of 'moSL' class.")
  }
  spobj = wrap.sphere(newdata)
  myp   = base::ncol(object$parameters$location)-1
  
  # INTERMEDIATE VALUES
  old.mu    = object$parameters$location
  old.mobj  = wrap.sphere(old.mu)
  old.sigma = object$parameters$scale
  old.pi    = object$parameters$proportion
  old.d2med = moSL.d2medmat(spobj, old.mobj)
  
  # COMPUTE, EXTRACT, AND RETURN
  fin.eta   = moSL.eta(old.d2med, old.sigma, old.pi, myp)
  output    = spmix.getcluster(fin.eta)
  return(output)
}

#  S3 METHOD : DENSITY
#' @rdname moSL
#' @concept sphere
#' @export
density.moSL <- function(object, newdata){
  # PREPARE
  if (!inherits(object, "moSL")){
    stop("* density.moSL : input is not an object of 'moSL' class.")
  } 
  spobj = wrap.sphere(newdata)
  myp   = base::ncol(object$parameters$location)-1
  
  # INTERMEDIATE VALUES
  old.mu    = object$parameters$location
  old.mobj  = wrap.sphere(old.mu)
  old.sigma = object$parameters$scale
  old.pi    = object$parameters$proportion
  old.d2med = moSL.d2medmat(spobj, old.mobj)
  
  # COMPUTE AND RETURN
  evaldensity = density_mixture_SL(old.d2med, old.sigma, old.pi, myp)
  return(evaldensity)
}
#' @keywords internal
#' @noRd
density_mixture_SL <- function(d2med, sigmas, props, p){
  N = base::nrow(d2med)
  K = base::ncol(d2med)
  
  # normalizing constant
  eval_constant <- function(sigma){
    myfunc <- function(r){
      return(exp(-r/sigma)*(sin(r)^(p-1)))
    }
    t1   = 2*(pi^(p/2))/(gamma(p/2))
    t2   = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
    return(t1*t2) 
  }
  
  vecCsig = rep(0,K)
  for (k in 1:K){
    vecCsig[k] = eval_constant(sigmas[k])
  }
  
  # evaluate the density
  mixdensity = array(0,c(N,K))
  for (n in 1:N){
    for (k in 1:K){
      mixdensity[n,k] = props[k]*exp(-d2med[n,k]/sigmas[k])/vecCsig[k]
    }
  }
  
  # evaluate the output
  output = as.vector(base::rowSums(mixdensity))
  return(output)
}