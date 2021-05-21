#' Spherical Normal Distribution
#' 
#' We provide tools for an isotropic spherical normal (SN) distributions on 
#' a \eqn{(p-1)}-sphere in \eqn{\mathbf{R}^p} for sampling, density evaluation, and maximum likelihood estimation 
#' of the parameters where the density is defined as
#' \deqn{f_SN(x; \mu, \lambda) = \frac{1}{Z(\lambda)} \exp \left( -\frac{\lambda}{2} d^2(x,\mu) \right)}
#' for location and concentration parameters \eqn{\mu} and \eqn{\lambda} respectively and the normalizing constant \eqn{Z(\lambda)}.
#' 
#' 
#' @param data data vectors in form of either an \eqn{(n\times p)} matrix or a length-\eqn{n} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param mu a length-\eqn{p} unit-norm vector of location.
#' @param log a logical; \code{TRUE} to return log-density, \code{FALSE} for densities without logarithm applied.
#' @param lambda a concentration parameter that is positive.
#' @param n the number of samples to be generated.
#' @param method an algorithm name for concentration parameter estimation.i
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-5).}
#' }
#' 
#' @return 
#' \code{dspnorm} gives a vector of evaluated densities given samples. \code{rspnorm} generates 
#' unit-norm vectors in \eqn{\mathbf{R}^p} wrapped in a list. \code{mle.spnorm} computes MLEs and returns a list 
#' containing estimates of location (\code{mu}) and concentration (\code{lambda}) parameters.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------------
#' #          Example with Spherical Normal Distribution
#' #
#' # Given a fixed set of parameters, generate samples and acquire MLEs.
#' # Especially, we will see the evolution of estimation accuracy.
#' # -------------------------------------------------------------------
#' ## DEFAULT PARAMETERS
#' true.mu  = c(1,0,0,0,0)
#' true.lbd = 5
#' 
#' ## GENERATE DATA N=1000
#' big.data = rspnorm(1000, true.mu, true.lbd)
#' 
#' ## ITERATE FROM 50 TO 1000 by 10
#' idseq = seq(from=50, to=1000, by=10)
#' nseq  = length(idseq)
#' 
#' hist.mu  = rep(0, nseq)
#' hist.lbd = rep(0, nseq)
#' 
#' for (i in 1:nseq){
#'   small.data = big.data[1:idseq[i]]          # data subsetting
#'   small.MLE  = mle.spnorm(small.data) # compute MLE
#'   
#'   hist.mu[i]  = acos(sum(small.MLE$mu*true.mu)) # difference in mu
#'   hist.lbd[i] = small.MLE$lambda
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(idseq, hist.mu,  "b", pch=19, cex=0.5, main="difference in location")
#' plot(idseq, hist.lbd, "b", pch=19, cex=0.5, main="concentration param")
#' abline(h=true.lbd, lwd=2, col="red")
#' par(opar)
#' }
#' 
#' 
#' @name spnorm
#' @concept distribution
#' @rdname spnorm
NULL

#' @rdname spnorm
#' @export
dspnorm <- function(data, mu, lambda, log=FALSE){
  ## PREPROCESSING
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  FNAME  = "dspnorm"
  mu     = check_unitvec(mu, FNAME)
  lambda = check_num_nonneg(lambda, FNAME)
  p      = length(mu) # dimension

  ## EVALUATION
  #   1. normalizing constant
  nconstant = dspnorm.constant(lambda, p)
  #   2. case branching
  if (is.vector(x)){
    logmux = auxsphere_log(mu, x)
    output = exp(-(lambda/2)*sum(logmux*logmux))
  } else {
    dvec   = as.vector(cppdist_int_1toN(mu, x));
    output = exp((-lambda/2)*(dvec^2))
    
    # nx     = nrow(x)
    # output = rep(0,nx)
    # for (i in 1:nx){
    #   logmux = aux_log(mu, as.vector(x[i,]))
    #   output[i] = exp(-(lambda/2)*sum(logmux*logmux))
    # }
  }
  #   3. scale by normalizing constant and RETURN
  if (log){
    return(log(output)-log(nconstant))
  } else {
    return(exp(log(output)-log(nconstant)))
  } 
}
#' @keywords internal
#' @noRd
dspnorm.spobj <- function(spobj, mu, lambda, log=FALSE){
  ## PREPROCESSING
  x      = sp2mat(spobj)
  FNAME  = "dspnorm"
  mu     = check_unitvec(mu, FNAME)
  lambda = check_num_nonneg(lambda, FNAME)
  p      = length(mu) # dimension
  
  ## EVALUATION
  #   1. normalizing constant
  nconstant = dspnorm.constant(lambda, p)
  #   2. case branching
  if (is.vector(x)){
    logmux = auxsphere_log(mu, x)
    output = exp(-(lambda/2)*sum(logmux*logmux))
  } else {
    dvec   = as.vector(cppdist_int_1toN(mu, x));
    output = exp((-lambda/2)*(dvec^2))
    
    # nx     = nrow(x)
    # output = rep(0,nx)
    # for (i in 1:nx){
    #   logmux = aux_log(mu, as.vector(x[i,]))
    #   output[i] = exp(-(lambda/2)*sum(logmux*logmux))
    # }
  }
  #   3. scale by normalizing constant and RETURN
  if (log){
    return(log(output)-log(nconstant))
  } else {
    return(exp(log(output)-log(nconstant)))
  } 
}

#' @rdname spnorm
#' @export
rspnorm <- function(n, mu, lambda){
  ## PREPROCESSING
  FNAME = "rspnorm"
  n      = max(1, round(n))
  mu     = check_unitvec(mu, FNAME)
  lambda = check_num_nonneg(lambda, FNAME)
  D      = length(mu) # dimension
  
  ## ITERATE or RANDOM
  if (lambda==0){
    output = array(0,c(n,D))
    for (i in 1:n){
      tgt = stats::rnorm(D)
      output[i,] = tgt/ sqrt(sum(tgt^2))
    }
  } else {
    #   1. tangent vectors
    vectors = array(0,c(n,D))
    for (i in 1:n){
      vectors[i,] = rspnorm.single(mu, lambda, D)
    }
    #   2. exponential mapping
    output = array(0,c(n,D))
    for (i in 1:n){
      output[i,] = auxsphere_exp(mu, as.vector(vectors[i,]))
    }  
  }
  
  ## RETURN
  samples = list()
  if (n==1){
    samples[[1]] = as.vector(output)
  } else {
    for (i in 1:n){
      samples[[i]] = as.vector(output[i,])
    }
  }
  return(samples)
}

#' @rdname spnorm
#' @export
mle.spnorm <- function(data, method=c("Newton","Halley","Optimize","DE"), ...){
  ## PREPROCESSING
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  pars   = list(...)
  pnames = names(pars)
  myiter = max(50, ifelse(("maxiter"%in%pnames), pars$maxiter, 50))
  myeps  = min(1e-6, max(0, ifelse(("eps"%in%pnames), as.double(pars$eps), 1e-6)))
  alldip = c("newton","halley","optimize","de")
  myway  = match.arg(tolower(method), alldip)
  
  ## STEP 1. INTRINSIC MEAN
  N = length(spobj$data)
  myweight = rep(1/N, N)
  opt.mean = as.vector(inference_mean_intrinsic(spobj$name, spobj$data, myweight, myiter, myeps)$mean)
  
  ## STEP 2. OPTIMAL LAMBDA
  opt.lambda = switch(myway,
                      "de"       = lambda_method_DE(x, opt.mean, myiter, myeps),
                      "optimize" = lambda_method_opt(x, opt.mean, myiter, myeps),
                      "newton"   = lambda_method_newton(x, opt.mean, myiter, myeps),
                      "halley"   = lambda_method_halley(x, opt.mean, myiter, myeps))
  
  ## RETURN
  output = list(mu=opt.mean, lambda=opt.lambda)
  return(output)
}




# all others --------------------------------------------------------------
#' @keywords internal
#' @noRd
lambda_method_halley <- function(data, mean, myiter, myeps){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  d1N = as.vector(auxsphere_dist_1toN(mean, data))
  C   = base::sum(d1N^2)
  
  # 3. compute a negative log-likelihood function to be minimized
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = (lambda*C)/2
    term2 = n*log(dspnorm.constant(lambda,D))
    return(term1+term2)
  }
  
  t0 = n*log(2*(pi^((D-1)/2))/gamma((D-1)/2))
  opt.fun.red <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      return(stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = (lambda*C)/2
    term2 = n*log(dspnorm.constant(lambda,D)) + t0
    return(term1+term2)
  }
  
  # 4. run Halley's Method
  # 4-1. try several inputs and rough start over a grid
  ntest = 25
  grid.lambda = base::exp(seq(from=-2,to=2,length.out=ntest))*stats::var(d1N)
  grid.values = rep(0,ntest)
  for (i in 1:ntest){
    grid.values[i] = opt.fun.red(grid.lambda[i])
  }
  xold = grid.lambda[which.min(grid.values)]
  
  # 4-2. run iterations
  maxiter = myiter
  abstol  = myeps
  for (i in 1:maxiter){
    
    # h = min(abs(xold), sqrteps)/8
    h = min(abs(xold)/2, 1e-4)
    
    # 5-point grid evaluation
    grr = opt.fun.red(xold+2*h)
    gr  = opt.fun.red(xold+h)
    g   = opt.fun.red(xold)
    gl  = opt.fun.red(xold-h)
    gll = opt.fun.red(xold-2*h)
    
    gd1 = (gr-gl)/(2*h)                       # first derivative
    gd2 = (gr-(2*g)+gl)/(h^2)                 # second derivative
    gd3 = (grr - 2*gr + 2*gl - gll)/(2*(h^3)) # third derivative
    
    xnew    = xold - ((2*gd1*gd2)/(2*(gd2^2) - gd1*gd3))
    # print(paste("iteration for Halley's : ",i," done with lambda=",xnew, sep=""))
    xinc    = abs(xnew-xold)
    xold    = xnew
    
    if (xinc < abstol){
      break
    }
  }
  
  # 5. return an optimal solution
  return(xold)
}
#' @keywords internal
#' @noRd
lambda_method_newton <- function(data, mean, myiter, myeps){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  d1N = as.vector(auxsphere_dist_1toN(mean, data))
  C   = base::sum(d1N^2)
  
  # 3. compute a negative log-likelihood function to be minimized
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = (lambda*C)/2
    term2 = n*log(dspnorm.constant(lambda,D))
    return(term1+term2)
  }
  
  t0 = n*log(2*(pi^((D-1)/2))/gamma((D-1)/2))
  opt.fun.red <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      return(stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = (lambda*C)/2
    term2 = n*log(dspnorm.constant(lambda,D)) + t0
    return(term1+term2)
  }
  
  # 4. run Newton's iteration
  # 4-1. try several inputs and rough start over a grid
  ntest = 25
  grid.lambda = base::exp(seq(from=-2,to=2,length.out=ntest))*stats::var(d1N)
  grid.values = rep(0,ntest)
  for (i in 1:ntest){
    grid.values[i] = opt.fun.red(grid.lambda[i])
  }
  xold = grid.lambda[which.min(grid.values)]
  
  # 4-2. run iterations
  maxiter = myiter
  for (i in 1:maxiter){
    # print(paste("iteration for Method 3 : ",i," initiated..", sep=""))
    h = min(abs(xold)/2, 1e-4)
    g.right = opt.fun.red(xold+h)
    g.mid   = opt.fun.red(xold)
    g.left  = opt.fun.red(xold-h)
    xnew    = xold - (h/2)*(g.right-g.left)/(g.right-(2*g.mid)+g.left)
    xinc    = abs(xnew-xold)
    xold    = xnew
    
    if (xinc < myeps){
      break
    }
  }
  
  # 5. return an optimal solution
  return(xold)
}
#' @keywords internal
#' @noRd
lambda_method_opt <- function(data, mean, myiter, myeps){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  d1N = as.vector(auxsphere_dist_1toN(mean, data))
  C   = base::sum(d1N^2)
  
  # 3. compute a log-likelihood function
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = -(lambda*C)/2
    term2 = -n*log(dspnorm.constant(lambda,D))
    return(term1+term2)
  }
  
  # 4. optimize a log-likelihood function with DEoptim
  myint  = c(0.01, 100)*stats::var(d1N)
  output = stats::optimize(opt.fun, interval=myint, maximum=TRUE, tol=myeps)$maximum
  return(output)
}
#' @keywords internal
#' @noRd
lambda_method_DE <- function(data, mean, myiter, myeps){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  d1N = as.vector(auxsphere_dist_1toN(mean, data))
  C   = base::sum(d1N^2)
  
  # 3. compute a of log-likelihood function
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = -(lambda*C)/2
    term2 = -n*log(dspnorm.constant(lambda,D))
    return(-(term1+term2))
  }
  
  # 4. optimize a log-likelihood function with DEoptim; careful with negative sign; take the last one as the optimum
  mymin  = stats::var(d1N)*0.01
  mymax  = stats::var(d1N)*100
  # output = as.double(tail(DEoptim(opt.fun, 10*.Machine$double.eps, 12345678, control=DEoptim.control(trace=FALSE))$member$bestmemit, n=1L))
  output = as.double(utils::tail(DEoptim::DEoptim(opt.fun, mymin, mymax, control=DEoptim.control(trace=FALSE, itermax=myiter, reltol=myeps))$member$bestmemit, n=1L))
  return(output)
}

#' @keywords internal
#' @noRd
sp2mat <- function(riemobj){
  N = length(riemobj$data)
  p = length(as.vector(riemobj$data[[1]]))
  
  output = array(0,c(N,p))
  for (n in 1:N){
    tmp = as.vector(riemobj$data[[n]])
    output[n,] = tmp/sqrt(sum(tmp^2))
  }
  return(output)
}
#' @keywords internal
#' @noRd
rspnorm.single <- function(mu, lambda, D){
  status = FALSE
  sqrts  = sqrt(1/(lambda + ((D-2)/pi))) # from the box : This is the correct one.
  #sqrts  = sqrt((lambda + ((D-2)/pi)))  # from the text
  while (status==FALSE){
    v = stats::rnorm(D,mean = 0, sd=sqrts)
    v = v-mu*sum(mu*v)       ## this part is something missed from the paper.
    v.norm = sqrt(sum(v^2))
    if (v.norm <= pi){
      if (v.norm <= sqrt(.Machine$double.eps)){
        r1 = exp(-(lambda/2)*(v.norm^2))  
      } else {
        r1 = exp(-(lambda/2)*(v.norm^2))*((sin(v.norm)/v.norm)^(D-2))
      }
      r2 = exp(-((v.norm^2)/2)*(lambda+((D-2)/pi)))
      r  = r1/r2
      u  = stats::runif(1)
      if (u <= r){
        status = TRUE
      }
    }
  }
  return(v)
}
#' @keywords internal
#' @noRd
auxsphere_log <- function(mu, x){
  # theta = base::acos(sum(x*mu))
  theta = tryCatch({base::acos(sum(x*mu))},
                   warning=function(w){
                     0
                   },error=function(e){
                     0
                   })
  if (abs(theta)<10*(.Machine$double.eps)){
    output = x-mu*(sum(x*mu))
  } else {
    output = (x-mu*(sum(x*mu)))*theta/sin(theta)
  }
  return(output)
}
#' @keywords internal
#' @noRd
auxsphere_exp <- function(x, d){
  nrm_td = norm(matrix(d),"f")
  if (nrm_td < sqrt(.Machine$double.eps)){
    output = x;
  } else {
    output = cos(nrm_td)*x + (sin(nrm_td)/nrm_td)*d; 
  }
  return(output)
}
#' @keywords internal
#' @noRd
auxsphere_dist_1toN <- function(x, maty){
  dist_one <- function(y){
    logxy = auxsphere_log(x, y)
    return(sqrt(sum((logxy)^2)))
  }
  return(as.vector(apply(maty, 1, dist_one)))
}
#' @keywords internal
#' @noRd
dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
  myfunc <- function(r){
    return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
  }
  t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
  t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
  return(t1*t2)
}



# # (EX1) COMPARE ALL METHODS
# myp   = 5
# mylbd = stats::runif(1, min=0.001, max=5)
# myn   = 2000
# mymu  = rnorm(myp)
# mymu  = mymu/sqrt(sum(mymu^2))
# myx   = Riemann::rspnorm(myn, mymu, lambda=mylbd)
# mle.spnorm(myx, method="de")
# mle.spnorm(myx, method="Optimize")
# mle.spnorm(myx, method="newton")
# mle.spnorm(myx, method="halley")
# mymu
# mylbd
# 
# # (EX2) COMPARE RUN TIME
# library(ggplot2)
# library(microbenchmark)  # time comparison of multiple methods
# lbdtime <- microbenchmark(
#   newton = mle.spnorm(myx, method="newton"),
#   halley = mle.spnorm(myx, method="halley"),
#   Roptim = mle.spnorm(myx, method="optimize"),
#   DEopt  = mle.spnorm(myx, method="DE"), times=10L
# )
# autoplot(lbdtime)