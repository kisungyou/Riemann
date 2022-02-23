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
    new.label = apply(new.eta, 1, whic.max)
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
  
  if (!same.lambda){
    #     par.k = ((myp-1)*myk) + myk + (myk-1)  
    #   } else {
    #     par.k = ((myp-1)*myk) + 1   + (myk-1)  
    #   }
    #   AIC = -2*old.loglkd + 2*par.k
    #   BIC = -2*old.loglkd + par.k*log(myn)
    #   HQIC = -2*old.loglkd + 2*par.k*log(log(myn))
    #   AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
    #   
    #   infov = c(AIC, AICc, BIC, HQIC)
    #   names(infov) = c("AIC","AICc","BIC","HQIC")
    #   
    #   ## RETURN
    #   output = list()
    #   output$cluster  = spmix.getcluster(old.eta)
    #   output$loglkd   = old.loglkd  # max loglkd
    #   output$criteria = infov       # min AIC/AICc/BIC/HQIC
    #   output$parameters = list(proportion=old.pi, center=old.mu, concentration=old.lambda)
    #   output$membership = old.eta
    #   return(structure(output, class="mixspnorm"))  
  
}




   
#   ## INFORMATION CRITERION
#   if (!same.lambda){
#     par.k = ((myp-1)*myk) + myk + (myk-1)  
#   } else {
#     par.k = ((myp-1)*myk) + 1   + (myk-1)  
#   }
#   AIC = -2*old.loglkd + 2*par.k
#   BIC = -2*old.loglkd + par.k*log(myn)
#   HQIC = -2*old.loglkd + 2*par.k*log(log(myn))
#   AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
#   
#   infov = c(AIC, AICc, BIC, HQIC)
#   names(infov) = c("AIC","AICc","BIC","HQIC")
#   
#   ## RETURN
#   output = list()
#   output$cluster  = spmix.getcluster(old.eta)
#   output$loglkd   = old.loglkd  # max loglkd
#   output$criteria = infov       # min AIC/AICc/BIC/HQIC
#   output$parameters = list(proportion=old.pi, center=old.mu, concentration=old.lambda)
#   output$membership = old.eta
#   return(structure(output, class="mixspnorm"))
# }



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
  output = as.double(stats::optimize(opt.fun, interval=myint, maximum=FALSE, tol=myeps)$minimum)
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



# auxiliary renewed -------------------------------------------------------

#  7. spmix.getcluster : find the label of the data given the membership
#  8. spmix.hard & spmix.stochastic




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


