#' Spherical Laplace Distribution
#' 
#' This is a collection of tools for learning with spherical Laplace (SL) distribution 
#' on a \eqn{(p-1)}-dimensional sphere in \eqn{\mathbf{R}^p} including sampling, density evaluation, and 
#' maximum likelihood estimation of the parameters. The SL distribution is characterized by the following 
#' density function,
#' \deqn{f_{SL}(x; \mu, \sigma) = \frac{1}{C(\sigma)} \exp \left( -\frac{d(x,\mu)}{\sigma}  \right)}
#' for location and scale parameters \eqn{\mu} and \eqn{\sigma} respectively.
#' 
#' @param data data vectors in form of either an \eqn{(n\times p)} matrix or a length-\eqn{n} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param mu a length-\eqn{p} unit-norm vector of location.
#' @param sigma a scale parameter that is positive.
#' @param n the number of samples to be generated.
#' @param log a logical; \code{TRUE} to return log-density, \code{FALSE} for densities without logarithm applied.
#' @param method an algorithm name for concentration parameter estimation. It should be one of \code{"Newton"}, \code{"Optimize"}, and \code{"DE"} (case-sensitive).
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{maximum number of iterations to be run (default:50).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-6).}
#' \item{use.exact}{a logical to use exact (\code{TRUE}) or approximate (\code{FALSE}) updating rules (default: \code{FALSE}).}
#' }
#' 
#' @return 
#' \code{dsplaplace} gives a vector of evaluated densities given samples. \code{rsplaplace} generates 
#' unit-norm vectors in \eqn{\mathbf{R}^p} wrapped in a list. \code{mle.splaplace} computes MLEs and returns a list 
#' containing estimates of location (\code{mu}) and scale (\code{sigma}) parameters.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------------
#' #          Example with Spherical Laplace Distribution
#' #
#' # Given a fixed set of parameters, generate samples and acquire MLEs.
#' # Especially, we will see the evolution of estimation accuracy.
#' # -------------------------------------------------------------------
#' ## DEFAULT PARAMETERS
#' true.mu  = c(1,0,0,0,0)
#' true.sig = 1
#' 
#' ## GENERATE A RANDOM SAMPLE OF SIZE N=1000
#' big.data = rsplaplace(1000, true.mu, true.sig)
#' 
#' ## ITERATE FROM 50 TO 1000 by 10
#' idseq = seq(from=50, to=1000, by=10)
#' nseq  = length(idseq)
#' 
#' hist.mu  = rep(0, nseq)
#' hist.sig = rep(0, nseq)
#' 
#' for (i in 1:nseq){
#'   small.data = big.data[1:idseq[i]]             # data subsetting
#'   small.MLE  = mle.splaplace(small.data)        # compute MLE
#'   
#'   hist.mu[i]  = acos(sum(small.MLE$mu*true.mu)) # difference in mu
#'   hist.sig[i] = small.MLE$sigma
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(idseq, hist.mu,  "b", pch=19, cex=0.5, 
#'      main="difference in location", xlab="sample size")
#' plot(idseq, hist.sig, "b", pch=19, cex=0.5, 
#'      main="scale parameter", xlab="sample size")
#' abline(h=true.sig, lwd=2, col="red")
#' par(opar)
#' }
#' 
#' @name splaplace
#' @concept distribution
#' @rdname splaplace
NULL

#' @rdname splaplace
#' @export
dsplaplace <- function(data, mu, sigma, log=FALSE){
  ## PREPROCESSING
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  FNAME  = "dsplaplace"
  mu     = check_unitvec(mu, FNAME)
  sigma  = check_num_nonneg(sigma, FNAME)
  p      = length(mu)-1 # dimension along with paper's notation
  
  ## EVALUATION
  #   1. normalizing constant
  nconstant = dsplaplace.constant(sigma, p)
  #   2. case branching
  if (is.vector(x)){
    logmux = auxsphere_log(mu, x)
    output = exp(-sum(logmux*logmux)/sigma)
  } else {
    dvec   = as.vector(cppdist_int_1toN(mu, x));
    output = exp(-dvec/sigma)
  }
  #   3. scale by normalizing constant and RETURN
  if (log){
    return(log(output)-log(nconstant))
  } else {
    return(exp(log(output)-log(nconstant)))
  } 
}

#' @rdname splaplace
#' @export
rsplaplace <- function(n, mu, sigma){
  ## PREPROCESSING
  FNAME = "rsplaplace"
  n     = max(1, round(n))
  mu    = check_unitvec(mu, FNAME)
  sigma = check_num_nonneg(sigma, FNAME)
  D     = length(as.vector(mu))
  
  ## ITERATE or RANDOM
  output = array(0,c(n,D))
  if (10*sigma > .Machine$double.xmax){ # random
    for (i in 1:n){
      tgt = stats::rnorm(D)
      output[i,] = tgt/ sqrt(sum(tgt^2))
    }
  } else {
    for (i in 1:n){
      output[i,] = rsplaplace.single(mu, sigma)
    }
  }
  
  ## RETURN
  samples = vector("list", length=n)
  if (n==1){
    samples[[1]] = as.vector(output)
  } else {
    for (i in 1:n){
      samples[[i]] = as.vector(output[i,])
    }
  }
  return(samples)
}


#' @rdname splaplace
#' @export
mle.splaplace <- function(data, method=c("DE","Optimize","Newton"), ...){
  ## PREPROCESSING
  spobj  = wrap.sphere(data)
  x      = sp2mat(spobj)
  pars   = list(...)
  pnames = names(pars)
  
  if ("maxiter"%in%pnames){
    myiter = max(10, round(pars$maxiter))
  } else {
    myiter = 50
  }
  if ("eps"%in%pnames){
    myeps = min(1e-6, max(0, as.double(pars$eps)))
  } else {
    myeps = 1e-6
  }
  myway = tolower(match.arg(method))
  if ("use.exact" %in% pnames){
    use_exact = as.logical(pars$use.exact)
  } else {
    use_exact = FALSE
  }
  
  ## STEP 1. INTRINSIC MEDIAN
  N = length(spobj$data)
  myweight   = rep(1/N, N)
  opt.median = as.vector(inference_median_intrinsic(spobj$name, spobj$data, myweight, myiter, myeps)$median)
  
  ## STEP 2. OPTIMAL SIGMA
  opt.sigma = switch(myway,
                     "newton"      = sigma_method_newton(x, opt.median, myiter, myeps, use_exact),
                     "optimize"    = sigma_method_opt(x, opt.median, myiter, myeps),
                     "de"          = sigma_method_DE(x, opt.median, myiter, myeps))
  
  ## RETURN
  output = list(mu=opt.median, sigma=opt.sigma)
  return(output)
}



# estimation of the scale parameters --------------------------------------
#' @keywords internal
#' @noRd
sigma_method_DE <- function(data, median, myiter, myeps){
  # 1. parameters
  p = length(median)-1  # dimension S^p
  n = nrow(data)        # number of data 
  
  # 2. compute a constant
  d1N  = as.vector(auxsphere_dist_1toN(median, data))
  Chat = mean(d1N)
  
  # 3. the objective function
  opt.fun <- function(sigma){
    norm_constant <- function(par_p, par_sigma){
      # term : surface
      t1 = 2*(pi^(par_p/2))/gamma(par_p/2)
      # term : intergration
      myfunc <- function(par_r){
        return(exp(-par_r/par_sigma)*(sin(par_r)^(par_p-1)))
      }
      t2 = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      # return
      return(t1*t2)
    }
    term1 = Chat/sigma
    term2 = log(norm_constant(p, sigma))
    return(term1+term2)
  }
  
  # 4. optimize : DEoptim 
  mymin = stats::var(d1N)*0.01
  mymax = stats::var(d1N)*100
  output = as.double(utils::tail(DEoptim::DEoptim(opt.fun, mymin, mymax, control=DEoptim.control(trace=FALSE, itermax=myiter, reltol=myeps))$member$bestmemit, n=1L))
  
  return(output)
}
#' @keywords internal
#' @noRd
sigma_method_opt <- function(data, median, myiter, myeps){
  # 1. parameters
  p = length(median)-1  # dimension S^p
  n = nrow(data)  # number of data 
  
  # 2. compute a constant
  d1N  = as.vector(auxsphere_dist_1toN(median, data))
  Chat = mean(d1N)
  
  # 3. the objective function
  opt.fun <- function(sigma){
    norm_constant <- function(par_p, par_sigma){
      # term : surface
      t1 = 2*(pi^(par_p/2))/gamma(par_p/2)
      # term : intergration
      myfunc <- function(par_r){
        return(exp(-par_r/par_sigma)*(sin(par_r)^(par_p-1)))
      }
      t2 = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      # return
      return(t1*t2)
    }
    term1 = Chat/sigma
    term2 = log(norm_constant(p, sigma))
    return(term1+term2)
  }
  
  # 4. optimize a log-likelihood function
  myint  = c(0.01, 100)*mean(d1N)
  output = stats::optimize(opt.fun, interval=myint, maximum=FALSE, tol=myeps)$minimum
  return(output)
}
#' @keywords internal
#' @noRd
sigma_method_newton <- function(data, median, myiter, myeps, myexact){
  if (myexact){
    return(sigma_method_newton_exact(data, median, myiter, myeps))
  } else {
    return(sigma_method_newton_approx(data, median, myiter, myeps))
  }
}
#' @keywords internal
#' @noRd
sigma_method_newton_exact <- function(data, median, myiter, myeps){
  # 1. parameters
  p = length(median)-1
  n = nrow(data)
  
  # 2. compute a constant
  d1N  = as.vector(auxsphere_dist_1toN(median, data))
  Chat = mean(d1N)
  
  # 3. integral evaluators
  integral_I0 <- function(sigma){
    # define an objective
    tgt_funI0 <- function(r){
      return(exp(-r/sigma)*(sin(r)^(p-1)))
    }
    # integrate
    output = as.double(stats::integrate(tgt_funI0, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol = sqrt(.Machine$double.eps))$value)
    return(output)
  }
  integral_I1 <- function(sigma){
    # define an objective
    tgt_funI1 <- function(r){
      t1 = (r/(sigma^2))
      t2 = exp(-r/sigma)*(sin(r)^(p-1))
      return(t1*t2)
    }
    # integrate
    output = as.double(stats::integrate(tgt_funI1, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol = sqrt(.Machine$double.eps))$value)
    return(output)
  }
  integral_I2 <- function(sigma){
    # define an objective
    tgt_funI2 <- function(r){
      t1 = ((r^2)/(sigma^4)) - ((2*r)/(sigma^3))
      t2 = exp(-r/sigma)*(sin(r)^(p-1))
      return(t1*t2)
    }
    # integrate
    output = as.double(stats::integrate(tgt_funI2, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol = sqrt(.Machine$double.eps))$value)
    return(output)
  }
  
  # 4. initialize
  ntest = 20
  grid.sigma = base::exp(seq(from=-1,to=1,length.out=ntest))*sqrt(stats::var(d1N))
  grid.value = rep(0,ntest)
  for (i in 1:ntest){
    sigma_now     = grid.sigma[i]
    grid.value[i] = (Chat/sigma_now) + log(integral_I0(sigma_now))
  }
  sigma_old = grid.sigma[which.min(grid.value)]
  sigma_new = 0
  
  # 5. Newton-Raphson update
  for (it in 1:myiter){
    # compute : quantities
    valI0 = integral_I0(sigma_old)
    valI1 = integral_I1(sigma_old)
    valI2 = integral_I2(sigma_old)
    
    # compute : rationals
    term_top = -(Chat/(sigma_old^2)) + (valI1/valI0)
    term_bot = (2*Chat/(sigma_old^3)) + ((valI0*valI2 - (valI1^2))/(valI0^2))
  
    # update
    sigma_new = sigma_old - term_top/term_bot
    sigma_inc = abs(sigma_old - sigma_new)/abs(sigma_old)
    sigma_old = sigma_new
    if (sigma_inc < myeps){
      break
    }
  }
  
  # return
  return(sigma_old)
}


#' @keywords internal
#' @noRd
sigma_method_newton_approx <- function(data, median, myiter, myeps){
  # 1. parameters
  p = length(median)-1  # dimension S^p
  n = nrow(data)  # number of data 
  
  # 2. compute a constant
  d1N  = as.vector(auxsphere_dist_1toN(median, data))
  Chat = mean(d1N)
  
  # 3. the objective function
  # fun_g <- function(sigma){
  #   norm_constant <- function(par_p, par_sigma){
  #     # term : surface
  #     t1 = 2*(pi^(par_p/2))/gamma(par_p/2)
  #     # term : intergration
  #     myfunc <- function(par_r){
  #       return(exp(-par_r/par_sigma)*(sin(par_r)^(par_p-1)))
  #     }
  #     t2 = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
  #     # return
  #     return(t1*t2)
  #   }
  #   term1 = Chat/sigma
  #   term2 = log(norm_constant(p, sigma))
  #   return(term1+term2)
  # }
  vec_r = seq(from=0, to=pi, length.out=1000)
  inc_r = vec_r[2]-vec_r[1]
  fun_g <- function(sigma){
    # approximate integral
    vec_f = exp(-vec_r/sigma)*(sin(vec_r)^(p-1))
    
    term1 = Chat/sigma
    term2 = log((2*sum(vec_f)-(vec_f[1] + vec_f[length(vec_r)]))*inc_r/2)
    return(term1+term2)
  }
  
  # 4. Newton's Method
  # start with a grid
  ntest = 20
  grid.sigma = base::exp(seq(from=-1,to=1,length.out=ntest))*sqrt(stats::var(d1N))
  grid.value = rep(0,ntest)
  for (i in 1:ntest){
    grid.value[i] = fun_g(grid.sigma[i])
  }
  
  # iterate
  h_true    = 1e-4
  sigma_old = grid.sigma[which.min(grid.value)]
  sigma_new = 0
  for (it in 1:myiter){
    # compute : quantities
    h      = min(h_true, sigma_old/2)
    eval_l = fun_g(sigma_old-h)
    eval_r = fun_g(sigma_old+h)
    eval_m = fun_g(sigma_old)
    
    # compute : derivatives
    gderiv1 = (eval_r - eval_l)/(2*h)
    gderiv2 = (eval_r - 2*eval_m + eval_l)/(h^2)
    # gderiv1 = (fun_g(sigma_old + h) - fun_g(sigma_old - h))/(2*h) 
    # gderiv2 = (fun_g(sigma_old+h) - 2*fun_g(sigma_old) + fun_g(sigma_old-h))/(2*(h^2))
    # gderiv1 = -Chat/(sigma_old^2) + valQ/valP
    # gderiv2 = 2*Chat/(sigma_old^3) + ((valP*valR - (valQ^2))/(valP^2))
    
    # update
    sigma_new = sigma_old - gderiv1/gderiv2
    sigma_inc = abs(sigma_old - sigma_new)/abs(sigma_old)
    sigma_old = sigma_new
    if (sigma_inc < myeps){
      break
    }
  }
  
  # return
  return(sigma_old)
}


# auxiliary functions for 'splaplace' distributions -----------------------
#  Formula for C_p(sigma)
#' @keywords internal
#' @noRd
dsplaplace.constant <- function(sigma, p){
  # define a function
  myfunc <- function(r){
    return(exp(-r/sigma)*(sin(r)^(p-1)))
  }
  
  # compute
  t1 = 2*(pi^(p/2))/(gamma(p/2))
  t2 = stats::integrate(myfunc, lower=sqrt(.Machine$double.eps), upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
  return(t1*t2)
}

#  Rejection sampling for the splaplace distribution
#' @keywords internal
#' @noRd
rsplaplace.single <- function(mu, sigma){
  status  = TRUE
  counter = 0
  while (status){
    # draw a single sample from 'spnormal'
    y = Riemann::rspnorm(1, mu, 1/sigma)[[1]]
    r = rsplaplace.dist(y, as.vector(mu))
    thr = exp(((r^2)/(2*sigma)) - (r/sigma) - (pi*(pi-2)/(2*sigma)))
    
    #counter = counter + 1
    if (stats::runif(1) < thr){
      #print(paste0("sampling done with sigma=", round(sigma, 5), " at iteration=",counter))
      status=FALSE
    }
    if (counter == 1e+4){
      #print(paste0("sampling stops with sigma=", round(sigma, 5), " at iteration=",counter))
      status=FALSE
    }
  }
  return(y)
}
#' @keywords internal
#' @noRd
rsplaplace.dist <- function(x, y){
  if (sqrt(sum((x-y)^2)) < 100*.Machine$double.eps){
    return(0)
  } else {
    return(acos(sum(x*y)))
  }
}

# ## IN-CODE TEST
# true.mu  = c(1,0,0,0,0)
# true.lbd = 0.01
# 
# ## GENERATE DATA N=1000
# small.data = rspnorm(1000, true.mu, true.lbd)
# 
# ## COMPARE FOUR METHODS
# test1 = mle.splaplace(small.data, method="Optimize")
# test2 = mle.splaplace(small.data, method="DE")
# test3 = mle.splaplace(small.data, method="Newton", use.exact=FALSE)
# test4 = mle.splaplace(small.data, method="Newton", use.exact=TRUE)


# 
# microbenchmark::microbenchmark(
#   test1 = mle.splaplace(small.data, method="optimize"),
#   test2 = mle.splaplace(small.data, method="DE"),
#   test3 = mle.splaplace(small.data, method="newton", eps=1e-4),
#   times = 5L)
