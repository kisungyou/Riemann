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
#' @param sigma a positive scale parameter; \code{Inf} gives the uniform limit.
#' @param n the number of samples to be generated.
#' @param log a logical; \code{TRUE} to return log-density, \code{FALSE} for densities without logarithm applied.
#' @param method an algorithm name for scale parameter estimation. It should be one of \code{"Newton"}, \code{"Optimize"}, and \code{"DE"} (case-sensitive).
#' @param ... extra parameters for computations, including\describe{
#' \item{maxiter}{iteration budget for each location, likelihood-search, or polishing stage; rounded and raised to at least 10 (default: 50).}
#' \item{eps}{positive tolerance, capped at 1e-6; the relative likelihood-score tolerance is additionally floored at 1e-10 (default: 1e-6).}
#' \item{use.exact}{for Newton, use moment-based derivatives (\code{TRUE}) or finite differences of the score (\code{FALSE}, default). Both use adaptive radial quadrature.}
#' }
#' 
#' @return 
#' \code{dsplaplace} gives a vector of evaluated densities given samples. \code{rsplaplace} generates 
#' unit-norm vectors in \eqn{\mathbf{R}^p} wrapped in a list. \code{mle.splaplace} computes MLEs and returns a list 
#' containing estimates of location (\code{mu}) and scale (\code{sigma}) parameters.
#' 
#' @details
#' Scale estimation uses an adaptively bracketed likelihood in the reciprocal
#' scale. Newton updates are safeguarded by that bracket; Optimize and DE search
#' the log reciprocal scale and use safeguarded Newton polishing if their score
#' has not reached the requested tolerance. The uniform boundary is returned as
#' \code{sigma = Inf}. For coincident observations the likelihood is unbounded
#' and \code{sigma = 0} is returned with a warning; this point-mass limit has no
#' surface-area density and is not accepted by \code{dsplaplace} or
#' \code{rsplaplace}. The location is obtained by local intrinsic-median
#' optimization and need not be globally optimal for data spread across the
#' sphere. Failure of the location iteration to converge produces a warning;
#' the returned scale then optimizes the likelihood conditional on that last
#' location estimate. Log densities use the log kernel and logarithmic
#' normalizer directly.
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
#' @section Validation status:
#' This retained legacy interface is experimental. Its full numerical and
#' statistical contract has not been independently verified across supported
#' inputs. See \code{\link{riem-method-contracts}} and the installed contract
#' table for method-specific assumptions, restrictions, and evidence scope.
dsplaplace <- function(data, mu, sigma, log=FALSE){
  x <- sp2mat(wrap.sphere(data))
  mu <- check_unitvec(mu, "dsplaplace")
  sigma <- splaplace_check_scale(sigma)
  dvec <- sphere_distribution_distances(mu, x)
  eta <- 1/sigma
  logdensity <- -dvec/sigma -
    sphere_radial_stats(eta, length(mu), 1, moments = FALSE)$logZ
  if (log) logdensity else exp(logdensity)
}

#' @keywords internal
#' @noRd
splaplace_check_scale <- function(sigma) {
  if (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0)
    stop("'sigma' must be positive (Inf denotes the uniform limit).", call. = FALSE)
  if (is.finite(sigma) && !is.finite(1/sigma))
    stop("'sigma' is too small to represent its reciprocal rate.", call. = FALSE)
  sigma
}

#' @rdname splaplace
#' @export
rsplaplace <- function(n, mu, sigma){
  ## PREPROCESSING
  FNAME = "rsplaplace"
  n     = max(1, round(n))
  mu    = check_unitvec(mu, FNAME)
  sigma = splaplace_check_scale(sigma)
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
  
  controls <- sphere_distribution_controls(pars)
  myiter <- controls$maxiter
  myeps <- controls$eps
  myway = tolower(match.arg(method))
  if ("use.exact" %in% pnames){
    if (!is.logical(pars$use.exact) || length(pars$use.exact) != 1L || is.na(pars$use.exact))
      stop("'use.exact' must be TRUE or FALSE.", call. = FALSE)
    use_exact = pars$use.exact
  } else {
    use_exact = FALSE
  }
  
  ## STEP 1. INTRINSIC MEDIAN
  opt.median <- sphere_distribution_location(spobj, x, myiter, myeps, median = TRUE)
  
  ## STEP 2. OPTIMAL SIGMA
  opt.sigma = switch(myway,
                     "newton"      = sigma_method_newton(x, opt.median, myiter, myeps, use_exact),
                     "optimize"    = sigma_method_opt(x, opt.median, myiter, myeps),
                     "de"          = sigma_method_DE(x, opt.median, myiter, myeps))
  
  ## RETURN
  output = list(mu=opt.median, sigma=opt.sigma)
  return(output)
}



# Scale estimation uses the natural rate eta = 1/sigma.
#' @keywords internal
#' @noRd
sigma_fit <- function(data, median, myiter, myeps, method, exact = TRUE) {
  d <- if (sphere_coincident_rows(data)) rep(0, nrow(data)) else
    sphere_distribution_distances(median, data)
  eta <- sphere_radial_mle(mean(d), length(median), 1, method,
                          myiter, myeps, exact)
  if (is.infinite(eta) && all(d == 0))
    warning("Coincident observations have no positive scale MLE; returning sigma = 0.", call. = FALSE)
  1/eta
}
#' @keywords internal
#' @noRd
sigma_method_DE <- function(data, median, myiter, myeps)
  sigma_fit(data, median, myiter, myeps, "de")
#' @keywords internal
#' @noRd
sigma_method_opt <- function(data, median, myiter, myeps)
  sigma_fit(data, median, myiter, myeps, "optimize")
#' @keywords internal
#' @noRd
sigma_method_newton <- function(data, median, myiter, myeps, myexact)
  sigma_fit(data, median, myiter, myeps, "newton", exact = myexact)
#' @keywords internal
#' @noRd
sigma_method_newton_exact <- function(data, median, myiter, myeps)
  sigma_fit(data, median, myiter, myeps, "newton", exact = TRUE)
#' @keywords internal
#' @noRd
sigma_method_newton_approx <- function(data, median, myiter, myeps)
  sigma_fit(data, median, myiter, myeps, "newton", exact = FALSE)
#' @keywords internal
#' @noRd
dsplaplace.constant <- function(sigma, p){
  exp(sphere_radial_stats(1/sigma, p+1, 1, moments = FALSE)$logZ)
}

#  Rejection sampling for the splaplace distribution
#' @keywords internal
#' @noRd
rsplaplace.single <- function(mu, sigma){
  # Try the rejection first
  tmp_reject = rsplaplace.single.rejection(mu, sigma)
  if (!tmp_reject$status){
    return(as.vector(tmp_reject$y))
  } else {
    return(rsplaplace.single.metropolis(mu, sigma))
  }
}

#' @keywords internal
#' @noRd
rsplaplace.single.rejection <- function(mu, sigma){
  status  = TRUE
  counter = 0
  while (status){
    # draw a single sample from 'spnormal'
    y = Riemann::rspnorm(1, mu, 1/sigma)[[1]]
    r = rsplaplace.dist(y, as.vector(mu))
    thr = exp(((r^2)/(2*sigma)) - (r/sigma) - (pi*(pi-2)/(2*sigma)))
    
    counter = counter + 1
    if (stats::runif(1) < thr){
      status=FALSE
    }
    if (counter >= 50){
      break
    }
  }
  
  output = list()
  output$y = y
  output$status = status 
  return(output)
}
#' @keywords internal
#' @noRd
rsplaplace.single.metropolis <- function(mu, sigma){
  myp  = length(mu)
  yold = mu + stats::rnorm(myp, mean=0, sd=sigma)
  yold = yold/sqrt(sum(yold^2))
  
  for (i in 1:50){
    # perturbation
    ytmp = yold + stats::rnorm(myp, mean=0, sd=sigma)
    ytmp = ytmp/sqrt(sum(ytmp^2))
    
    # compute a ratio
    threshold = exp((rsplaplace.dist(yold, mu)-rsplaplace.dist(ytmp, mu))/sigma)
    if (stats::runif(1) < threshold){
      ynew = ytmp
    } else {
      ynew = yold
    }
    yold = ynew
  }
  return(yold)
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
