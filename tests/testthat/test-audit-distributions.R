# Independent circle identities: these do not use package quadrature.
sn_circle_logZ <- function(lambda) {
  if (lambda == 0) return(log(2*pi))
  .5*log(2*pi/lambda) + log(2*stats::pnorm(pi*sqrt(lambda))-1)
}
sl_circle_logZ <- function(sigma) {
  if (is.infinite(sigma)) return(log(2*pi))
  log(2) + log(sigma) + log(-expm1(-pi/sigma))
}
circle_data <- function(angles) cbind(cos(angles), sin(angles))

test_that("spherical log densities preserve finite values below density range", {
  expect_equal(dspnorm(rbind(c(0,1)), c(1,0), 1000, log=TRUE),
               -1000*(pi/2)^2/2 - sn_circle_logZ(1000), tolerance=1e-12)
  expect_equal(dsplaplace(rbind(c(-1,0)), c(1,0), .002, log=TRUE),
               -pi/.002 - sl_circle_logZ(.002), tolerance=1e-12)
  expect_identical(dspnorm(rbind(c(0,1)), c(1,0), 1000), 0)
  expect_identical(dsplaplace(rbind(c(-1,0)), c(1,0), .002), 0)
  angles <- c(0, 1e-10, .1, 1, pi)
  for (lambda in c(0, 1e-4, 1, 1000, 1e100)) {
    expected <- -lambda*angles^2/2 - sn_circle_logZ(lambda)
    actual <- dspnorm(circle_data(angles), c(1,0), lambda, log=TRUE)
    expect_equal(actual, expected, tolerance=1e-10)
    expect_equal(Riemann:::dspnorm.spobj(wrap.sphere(circle_data(angles)),
                 c(1,0), lambda, log=TRUE), actual)
  }
  for (sigma in c(Inf, 1000, 1, .002, 1e-100)) {
    expected <- -angles/sigma - sl_circle_logZ(sigma)
    expect_equal(dsplaplace(circle_data(angles), c(1,0), sigma, log=TRUE),
                 expected, tolerance=1e-10)
  }
  nonaxis <- rep(1/sqrt(3),3)
  expect_equal(dspnorm(matrix(nonaxis,1),nonaxis,1e100,log=TRUE),
               -log(2*pi/1e100), tolerance=1e-10)
  expect_error(dsplaplace(rbind(c(1,0)), c(1,0), 0), "positive")
  expect_error(rsplaplace(1, c(1,0), 0), "positive")
  expect_error(dspnorm(rbind(c(1,0)), c(1,0), Inf), "nonnegative")
  expect_error(dspnorm(rbind(c(1,0)), c(1,0,0), 1), "same dimension")
})

test_that("radial normalization agrees with separate angular integration", {
  for (D in c(2,3,8)) {
    mu <- c(1, rep(0,D-1))
    surface <- 2*pi^((D-1)/2)/gamma((D-1)/2)
    for (lambda in c(0, .3, 4)) {
      reference <- surface * integrate(function(r)
        exp(-lambda*r*r/2)*sin(r)^(D-2), 0, pi, rel.tol=1e-11)$value
      expect_equal(dspnorm(matrix(mu,1),mu,lambda,log=TRUE),
                   -log(reference), tolerance=1e-9)
    }
    for (sigma in c(.2, 3, Inf)) {
      reference <- surface * integrate(function(r)
        exp(-r/sigma)*sin(r)^(D-2), 0, pi, rel.tol=1e-11)$value
      expect_equal(dsplaplace(matrix(mu,1),mu,sigma,log=TRUE),
                   -log(reference), tolerance=1e-9)
    }
  }
  expect_equal(integrate(function(a) dspnorm(circle_data(a),c(1,0),3),
                         -pi,pi,rel.tol=1e-9)$value, 1, tolerance=1e-8)
  expect_equal(integrate(function(a) dsplaplace(circle_data(a),c(1,0),.4),
                         -pi,pi,rel.tol=1e-9)$value, 1, tolerance=1e-8)
  # These high-concentration log normalizers remain representable even when
  # the surface normalizer itself underflows (nine intrinsic dimensions).
  D <- 10
  lambda <- 1e100
  normal.limit <- ((D-1)/2)*log(2*pi/lambda)
  laplace.limit <- log(2) + ((D-1)/2)*log(pi)-lgamma((D-1)/2) +
    lgamma(D-1) + (D-1)*log(1e-100)
  mu <- c(1,rep(0,D-1))
  expect_equal(dspnorm(matrix(mu,1),mu,lambda,log=TRUE),
               -normal.limit, tolerance=1e-10)
  expect_equal(dsplaplace(matrix(mu,1),mu,1e-100,log=TRUE),
               -laplace.limit, tolerance=1e-10)
})

test_that("all normal concentration methods reach the circle likelihood optimum", {
  for (a in list(c(-.2,-.15,-.1,-.05,.05,.1,.15,.2), c(-.1,.1),
                 1e-5*c(-.2,-.1,.1,.2))) {
    # At these concentrations the circle truncation tail is negligible.
    reference <- 1/mean(a^2)
    for (method in c("Newton","Halley","Optimize","DE")) {
      set.seed(104)
      fit <- mle.spnorm(circle_data(a), method=method, maxiter=100, eps=1e-8)
      expect_equal(fit$mu, c(1,0), tolerance=1e-8)
      expect_equal(fit$lambda/reference, 1, tolerance=2e-7)
    }
  }
})

test_that("all Laplace scale methods reach the circle likelihood optimum", {
  for (a in list(c(-.02,-.015,-.01,-.005,.005,.01,.015,.02),
                 c(-.1,.1), 1e-5*c(-.02,-.01,.01,.02))) {
    reference <- mean(abs(a))
    for (method in c("Newton","Optimize","DE")) {
      for (exact in if (method == "Newton") c(FALSE,TRUE) else TRUE) {
        set.seed(114)
        fit <- mle.splaplace(circle_data(a), method=method, maxiter=100,
                             eps=1e-8, use.exact=exact)
        # An even sample has a whole median arc between its middle angles.
        # A two-point sample may validly return either endpoint of that arc.
        fitted.angle <- atan2(fit$mu[2],fit$mu[1])
        middle <- sort(a)[c(length(a)/2,length(a)/2+1)]
        expect_gte(fitted.angle, middle[1]-1e-8)
        expect_lte(fitted.angle, middle[2]+1e-8)
        expect_equal(mean(abs(a-fitted.angle)),reference,tolerance=1e-8)
        expect_equal(fit$sigma/reference, 1, tolerance=2e-7)
      }
    }
  }
})

test_that("diffuse circle likelihoods retain their nonnegligible tail terms", {
  a <- c(-2.8,-1.4,1.4,2.8)
  for (method in c("Newton","Halley","Optimize","DE")) {
    set.seed(204)
    fit <- mle.spnorm(circle_data(a),method=method,maxiter=100,eps=1e-8)
    d <- abs(atan2(sin(a-atan2(fit$mu[2],fit$mu[1])),
                  cos(a-atan2(fit$mu[2],fit$mu[1]))))
    ref <- optimize(function(lambda) lambda*mean(d^2)/2+sn_circle_logZ(lambda),
                    c(1e-6,10),tol=1e-10)$minimum
    expect_equal(fit$lambda,ref,tolerance=2e-7)
  }
  for (method in c("Newton","Optimize","DE")) {
    set.seed(214)
    fit <- mle.splaplace(circle_data(a),method=method,maxiter=100,eps=1e-8)
    d <- abs(atan2(sin(a-atan2(fit$mu[2],fit$mu[1])),
                  cos(a-atan2(fit$mu[2],fit$mu[1]))))
    ref <- optimize(function(sigma) mean(d)/sigma+sl_circle_logZ(sigma),
                    c(.01,100),tol=1e-10)$minimum
    expect_equal(fit$sigma,ref,tolerance=2e-7)
  }
})

test_that("radial MLE boundaries are explicit for every estimation method", {
  # Conditional radial likelihood at a supplied location, deliberately not a
  # claim that this location is a global mean for diffuse observations.
  antipodes <- rbind(c(1,0),c(-1,0))
  for (method in c("newton","halley","optimize","de")) {
    expect_identical(Riemann:::lambda_fit(antipodes,c(1,0),50,1e-8,method),0)
  }
  for (method in c("newton","optimize","de")) {
    expect_identical(Riemann:::sigma_fit(antipodes,c(1,0),50,1e-8,method),Inf)
  }
  for (mu in list(c(1,0), rep(1/sqrt(3),3))) {
    data <- matrix(rep(mu,3),3,byrow=TRUE)
    for (method in c("Newton","Halley","Optimize","DE")) {
      expect_warning(fit <- mle.spnorm(data,method=method),"no finite concentration MLE")
      expect_identical(fit$lambda,Inf)
      expect_equal(fit$mu,mu,tolerance=1e-14)
    }
    for (method in c("Newton","Optimize","DE")) {
      expect_warning(fit <- mle.splaplace(data,method=method),"no positive scale MLE")
      expect_identical(fit$sigma,0)
      expect_equal(fit$mu,mu,tolerance=1e-14)
    }
  }
  expect_error(mle.spnorm(circle_data(c(-.1,.1)),eps=0),"positive finite")
  expect_error(mle.splaplace(circle_data(c(-.1,.1)),maxiter=Inf),"positive finite")
})


test_that("spherical radial scores recover independently specified rates", {
  # Set an exact population moment by independent unscaled angular quadrature,
  # then invert it. This exercises curvature factors beyond the circle.
  for (D in c(3,8)) for (power in c(1,2)) for (eta in c(.3,5)) {
    kernel <- function(r) exp(-eta*r^power)*sin(r)^(D-2)
    norm <- integrate(kernel,0,pi,rel.tol=1e-11)$value
    target <- integrate(function(r) r^power*kernel(r),0,pi,
                        rel.tol=1e-11)$value/norm
    for (method in c("newton","halley","optimize")) {
      estimate <- Riemann:::sphere_radial_mle(target,D,power,method,100,1e-8)
      expect_equal(estimate,eta,tolerance=2e-7)
    }
  }
  # On S^2, the Laplace radial normalizer has this independent closed form.
  for (eta in c(0,.3,5,1000)) {
    reference <- log(2*pi)+log1p(exp(-pi*eta))-log1p(eta^2)
    expect_equal(Riemann:::sphere_radial_stats(eta,3,1,FALSE)$logZ,
                 reference,tolerance=1e-10)
  }
})


test_that("Laplace likelihood rates use the full representable positive range", {
  # On the circle this is the exponential mean, with an unrepresentably small
  # truncation correction. The rate exceeds double.xmax/2 but remains finite.
  rate <- Riemann:::sphere_radial_mle(1e-308,2,1,"newton",100,1e-8)
  expect_true(is.finite(rate))
  expect_equal(rate/1e308,1,tolerance=1e-8)
})


test_that("the finite normal concentration boundary does not overflow", {
  eta <- Riemann:::sphere_radial_mle(1/.Machine$double.xmax,2,2,
                                    "newton",100,1e-8)
  expect_true(is.finite(2*eta))
  expect_equal(2*eta/.Machine$double.xmax,1,tolerance=1e-8)
})

test_that("spherical fits report an unconverged location solve", {
  set.seed(1)
  x <- matrix(rnorm(100),20,5)
  x <- x/sqrt(rowSums(x^2))
  for (fun in list(mle.spnorm,mle.splaplace)) {
    expect_warning(fit <- fun(x,method="Newton",maxiter=10),
                   "location iteration did not converge")
    expect_true(all(is.finite(fit$mu)))
    expect_equal(sum(fit$mu^2),1,tolerance=1e-12)
  }
})
