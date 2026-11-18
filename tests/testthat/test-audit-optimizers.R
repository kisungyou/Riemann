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

test_that("MACG search samples remain valid as elite subspaces concentrate", {
  # The exact trace maximum on Gr(2,4) is the sum of the top two eigenvalues.
  A <- diag(c(5, 3, 1, -2))
  objective <- function(q) -sum(q * (A %*% q))
  set.seed(17)
  fit <- grassmann.optmacg(objective, 4, 2, n.start = 2, popsize = 48,
                           ratio = .25, maxiter = 30)
  expect_equal(crossprod(fit$solution), diag(2), tolerance = 1e-12)
  expect_equal(fit$cost, objective(fit$solution))
  expect_gte(fit$cost, -8 - 1e-12)
  expect_lt(fit$cost, -7.9)

  # One elite frame has a singular span for k < p. Its sampling proposal
  # must remain positive definite without modifying the public MACG estimator.
  calls <- 0L
  observed <- numeric()
  set.seed(27)
  single <- grassmann.optmacg(function(q) {
    calls <<- calls + 1L
    cost <- objective(q)
    observed <<- c(observed, cost)
    cost
  }, 4, 2, n.start = 2, popsize = 4, ratio = .1, maxiter = 2)
  expect_equal(calls, 24L) # two generations per run plus an eight-frame pilot
  expect_identical(single$cost, min(observed))
  expect_equal(crossprod(single$solution), diag(2), tolerance = 1e-12)
  expect_equal(single$cost, objective(single$solution))
})

test_that("the full Grassmann subspace has a single attainable objective", {
  evaluations <- 0L
  fit <- grassmann.optmacg(function(x) {
    evaluations <<- evaluations + 1L
    sum(x^2)
  }, p = 3, k = 3)
  expect_identical(evaluations, 1L)
  expect_equal(fit$solution, diag(3))
  expect_equal(fit$cost, 3)
})
