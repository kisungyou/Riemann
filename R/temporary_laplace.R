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
  
  
  old.lambda = spmix.updatelbd(old.eta, old.d2mat, myp, homogeneous = same.lambda)
  old.pi     = as.vector(base::colSums(old.eta)/myn)
  old.loglkd = spmix.loglkd(old.d2mat, old.lambda, old.pi, myp)
}


mixspnorm <- function(data, k=2, same.lambda=FALSE, variants=c("soft","hard","stochastic"), ...){
  ## PREPROCESSING
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
    
    # we need to stop if there is empty cluster
    new.label = apply(new.eta, 1, which.max)
    if (length(unique(new.label)) < myk){
      break
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



# auxiliary functions -----------------------------------------------------
# 1. moSL.median   : weighted Frechet medians
# 2. moSL.d2medmat : compute data-median pairwise distance

# 1. compute weighted Frechet median : deals with Riemdata object
#' @keywords internal
#' @noRd
moSL.median <- function(spobj, membership){
  N = base::length(spobj$data)
  P = base::length(as.vector(spobj$data[[1]]))
  K = base::ncol(membership)
  
  output = array(0,c(K,P))
  for (k in 1:K){
    parweight = as.vector(membership[,k])
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


# auxiliary renewed -------------------------------------------------------
#  3. spmix.solvelbd   : solve min A*lbd + B*log(Z(lbd))
#  4. spmix.updatelbd  : update lambda parameters
#  5. spmix.loglkd     : compute log-likelihood
#  6. spmix.eta        : compute eta (membership)
#  7. spmix.getcluster : find the label of the data given the membership
#  8. spmix.hard & spmix.stochastic




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
  
  vecZlbd = rep(0,K)
  for (k in 1:K){
    vecZlbd[k] = evalZ(lambdas[k])
  }
  # logprop = base::log(props)
  # for (k in 1:K){
  #   logZlbd[k] = base::log(evalZ(lambdas[k]))
  # }
  
  # evaluate the density
  mixdensity = array(0,c(N,K))
  for (n in 1:N){
    for (k in 1:K){
      mixdensity[n,k] = props[k]*exp(-(d2mat[n,k])*lambdas[k]/2)/vecZlbd[k]
    }
  }
  
  # evaluate the output
  output = base::sum(base::log(base::rowSums(mixdensity)))
  return(output)
  # term1 = base::sum(logprop)*N
  # term2 = -base::sum(logZlbd)*N
  # term3 = base::sum(d2mat%*%(-base::diag(lambdas)/2))
  # return(term1+term2+term3)
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


# xx = c(rspnorm(50, c(1,0,0), 20),
#        rspnorm(50, c(-1/sqrt(2),1/sqrt(2),0), 20),
#        rspnorm(50, c(0,0,-1), 20))
# 
# cc = rbind(mixspnorm(xx, k=2)$criteria,
# mixspnorm(xx, k=3)$criteria,
# mixspnorm(xx, k=4)$criteria,
# mixspnorm(xx, k=5)$criteria,
# mixspnorm(xx, k=6)$criteria,
# mixspnorm(xx, k=7)$criteria)
# 
# matplot(2:7, cc, type="b")


