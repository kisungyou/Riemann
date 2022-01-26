#' Spherical Laplace Distribution
#' 
#' @param mu a length-\eqn{p} unit-norm vector of location.
#' @param sigma a scale parameter that is positive.
#' @param n the number of samples to be generated.
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
mle.splaplace <- function(data, method=c("Newton","Halley","Optimize","DE"), ...){
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
  t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
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
    return(base::acos(base::sum(x*y)))
  }
}