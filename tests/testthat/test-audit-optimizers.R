test_that("Stiefel QR projection preserves frame orientation and is local", {
  for (frame in list(matrix(c(1,0,0),3,1), matrix(c(-1,0,0),3,1),
                     diag(3),diag(c(-1,1,1)))) {
    expect_equal(Riemann:::stiefel_qr_retract(frame),frame,tolerance=1e-14)
    set.seed(107)
    perturbation <- matrix(rnorm(length(frame),sd=1e-8),nrow(frame))
    projected <- Riemann:::stiefel_qr_retract(frame+perturbation)
    expect_lt(max(abs(projected-frame)),1e-6)
    expect_equal(crossprod(projected),diag(ncol(frame)),tolerance=1e-14)
  }
  # A Gaussian sign correction has full orientation support, including both
  # components of O(3); this checks support, not stochastic global convergence.
  set.seed(108)
  frames <- replicate(100,Riemann:::stiefel_qr_retract(matrix(rnorm(9),3)))
  signs <- vapply(seq_len(100),function(i) sign(det(frames[,,i])),0)
  expect_true(any(signs < 0))
  expect_true(any(signs > 0))
})

test_that("Stiefel public annealing evaluates both hemispheres", {
  seen <- numeric()
  objective <- function(x) { seen <<- c(seen,x[1,1]); -x[1,1] }
  set.seed(107)
  fit <- stiefel.optSA(objective,p=3,k=1,n.start=5,maxiter=100)
  expect_true(any(seen > 0))
  expect_true(any(seen < 0))
  expect_equal(crossprod(fit$solution),matrix(1),tolerance=1e-12)
  expect_identical(fit$cost,-fit$solution[1,1])
  # Tiny proposals must remain near a supplied positive frame, rather than
  # being reflected immediately to the negative hemisphere by QR signs.
  seen <- numeric()
  set.seed(109)
  fit <- stiefel.optSA(objective,p=3,k=1,init.val=matrix(c(1,0,0),3,1),
                       n.start=1,maxiter=10,stepsize=1e-8)
  expect_true(all(seen > .999))
  expect_error(stiefel.optSA(function(x) NA_real_,3,1),"finite numeric scalar")
})

test_that("Grassmann optimizer returns attained costs independent of offsets", {
  run <- function(offset) {
    set.seed(106)
    grassmann.optmacg(function(x) offset,p=2,k=1,n.start=10,
                      maxiter=20,popsize=50,ratio=.5)
  }
  low <- run(2)
  high <- run(2e7)
  negative <- run(-2e7)
  for (fit in list(low,high,negative)) {
    expect_true(is.matrix(fit$solution))
    expect_identical(dim(fit$solution),c(2L,1L))
    expect_equal(crossprod(fit$solution),matrix(1),tolerance=1e-12)
  }
  expect_identical(low$cost,2)
  expect_identical(high$cost,2e7)
  expect_identical(negative$cost,-2e7)
  expect_identical(low$solution,high$solution)
  expect_identical(low$solution,negative$solution)
  expect_error(grassmann.optmacg(function(x) Inf,2,1),"finite numeric scalar")
  expect_error(grassmann.optmacg(function(x) c(1,2),2,1),"finite numeric scalar")
})
