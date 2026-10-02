#' Spherical Normal Distribution
#' 
#' We provide tools for an isotropic spherical normal (SN) distributions on 
#' a \eqn{(p-1)}-sphere in \eqn{\mathbf{R}^p} for sampling, density evaluation, and maximum likelihood estimation 
#' of the parameters where the density is defined as
#' \deqn{f_{SN}(x; \mu, \lambda) = \frac{1}{Z(\lambda)} \exp \left( -\frac{\lambda}{2} d^2(x,\mu) \right)}
#' for location and concentration parameters \eqn{\mu} and \eqn{\lambda} respectively and the normalizing constant \eqn{Z(\lambda)}.
#' 
#' 
#' @param data data vectors in form of either an \eqn{(n\times p)} matrix or a length-\eqn{n} list.  See \code{\link{wrap.sphere}} for descriptions on supported input types.
#' @param mu a length-\eqn{p} unit-norm vector of location.
#' @param log a logical; \code{TRUE} to return log-density, \code{FALSE} for densities without logarithm applied.
#' @param lambda a finite nonnegative concentration parameter; zero is the uniform distribution.
#' @param n the number of samples to be generated.
#' @param method an algorithm name for concentration parameter estimation. It should be one of \code{"Newton"},\code{"Halley"},\code{"Optimize"}, and \code{"DE"} (case sensitive).
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{iteration budget for each location, likelihood-search, or polishing stage; rounded and raised to at least 10 (default: 50).}
#' \item{eps}{positive tolerance, capped at 1e-6; the relative likelihood-score tolerance is additionally floored at 1e-10 (default: 1e-6).}
#' }
#' 
#' @return 
#' \code{dspnorm} gives a vector of evaluated densities given samples. \code{rspnorm} generates 
#' unit-norm vectors in \eqn{\mathbf{R}^p} wrapped in a list. \code{mle.spnorm} computes MLEs and returns a list 
#' containing estimates of location (\code{mu}) and concentration (\code{lambda}) parameters.
#' 
#' @details
#' Concentration estimation uses an adaptively bracketed likelihood score.
#' Newton and Halley updates are safeguarded by that bracket; Optimize and DE
#' search the log concentration and use safeguarded Newton polishing if their
#' score has not reached the requested tolerance. All methods check the uniform
#' boundary, returning \code{lambda = 0} when appropriate. For coincident
#' observations the likelihood is unbounded and \code{lambda = Inf} is returned
#' with a warning; this point-mass limit has no density with respect to spherical
#' surface area and is not accepted by \code{dspnorm} or \code{rspnorm}.
#' The location is obtained by local intrinsic-mean optimization; a globally
#' optimal location is not guaranteed for data spread across the sphere.
#' Failure of the location iteration to converge produces a warning; the
#' returned concentration then optimizes the likelihood conditional on that
#' last location estimate.
#' Log densities are evaluated directly, including a logarithmic normalizer.
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
#' @references 
#' \insertRef{hauberg_2018_DirectionalStatisticsSpherical}{Riemann}
#' 
#' \insertRef{you_2022_ParameterEstimationModelbased}{Riemann}
#' 
#' @name spnorm
#' @concept distribution
#' @rdname spnorm
NULL

#' @rdname spnorm
#' @export
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
dspnorm <- function(data, mu, lambda, log=FALSE){
  dspnorm.spobj(wrap.sphere(data), mu, lambda, log)
}
#' @keywords internal
#' @noRd
dspnorm.spobj <- function(spobj, mu, lambda, log=FALSE){
  x <- sp2mat(spobj)
  mu <- check_unitvec(mu, "dspnorm")
  lambda <- check_num_nonneg(lambda, "dspnorm")
  dvec <- sphere_distribution_distances(mu, x)
  logdensity <- -(lambda/2)*dvec^2 -
    sphere_radial_stats(lambda/2, length(mu), 2, moments = FALSE)$logZ
  if (log) logdensity else exp(logdensity)
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
  
  controls <- sphere_distribution_controls(pars)
  myiter <- controls$maxiter
  myeps <- controls$eps
  myway = tolower(match.arg(method))
  
  ## STEP 1. INTRINSIC MEAN
  opt.mean <- sphere_distribution_location(spobj, x, myiter, myeps)
  
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




# Concentration methods share a bracketed likelihood and stable moments.
#' @keywords internal
#' @noRd
lambda_fit <- function(data, mean, myiter, myeps, method) {
  d <- if (sphere_coincident_rows(data)) rep(0, nrow(data)) else
    sphere_distribution_distances(mean, data)
  target <- base::mean(d^2)
  if (target == 0 && any(d > 0))
    warning("The finite concentration MLE exceeds the representable parameter range.", call. = FALSE)
  eta <- sphere_radial_mle(target, length(mean), 2, method, myiter, myeps)
  if (is.infinite(eta) && all(d == 0))
    warning("Coincident observations have no finite concentration MLE; returning lambda = Inf.", call. = FALSE)
  2*eta
}
#' @keywords internal
#' @noRd
lambda_method_halley <- function(data, mean, myiter, myeps)
  lambda_fit(data, mean, myiter, myeps, "halley")
#' @keywords internal
#' @noRd
lambda_method_newton <- function(data, mean, myiter, myeps)
  lambda_fit(data, mean, myiter, myeps, "newton")
#' @keywords internal
#' @noRd
lambda_method_opt <- function(data, mean, myiter, myeps)
  lambda_fit(data, mean, myiter, myeps, "optimize")
#' @keywords internal
#' @noRd
lambda_method_DE <- function(data, mean, myiter, myeps)
  lambda_fit(data, mean, myiter, myeps, "de")

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
  sphere_distribution_distances(x, maty)
}
#' @keywords internal
#' @noRd
dspnorm.constant <- function(lbd, D){
  exp(sphere_radial_stats(lbd/2, D, 2, moments = FALSE)$logZ)
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
