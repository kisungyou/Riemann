library(Riemann)

test_that("wrappers work", {
  
  ## 01sphere 
  data_sphere = array(0,c(5,3))
  for (i in 1:5){
    tgt = stats::rnorm(3)
    data_sphere[i,] = tgt/sqrt(sum(tgt^2))
  }
  riem01 = Riemann::wrap.sphere(data_sphere)
  
  ## 02spd
  data_spd = array(0,c(3,3,5))
  for (i in 1:5){
    dat = matrix(rnorm(10*3), ncol=3)
    data_spd[,,i] = stats::cov(dat)
  }
  riem02 = Riemann::wrap.spd(data_spd)
  
  ## 03correlation
  data_corr = array(0,c(3,3,5))
  for (i in 1:5){
    dat = matrix(rnorm(10*3), ncol=3)
    data_corr[,,i] = stats::cor(dat)
  }
  riem03 = Riemann::wrap.correlation(data_corr)
  
  ## 04stiefel & 05grassmann
  data_stiefel = array(0,c(4,2,5))
  for (i in 1:5){
    dat = matrix(rnorm(4*2), ncol=2)
    data_stiefel[,,i] = qr.Q(qr(dat))
  }
  riem04 = Riemann::wrap.stiefel(data_stiefel)
  riem05 = Riemann::wrap.grassmann(data_stiefel)
  
  ## 06rotation
  data_rotation = array(0,c(3,3,5))
  for (i in 1:5){
    data_rotation[,,i] = qr.Q(qr(matrix(rnorm(9), ncol=3)))
  }
  riem06 = Riemann::wrap.rotation(data_rotation)
  
  ## 07multinomial
  data_multinomial = array(0,c(5,3))
  for (i in 1:5){
    dat = abs(stats::rnorm(3))
    data_multinomial[i,] = dat/base::sum(dat)
  }
  riem07 = Riemann::wrap.multinomial(data_multinomial)
  
  ## 08spdk
  data_spdk = array(0,c(5,5,3))
  for (i in 1:3){
    data_spdk[,,i] = stats::cov(matrix(rnorm(10*5), ncol=5))
  }
  riem08 = Riemann::wrap.spdk(data_spdk, k=3)
  
  ## 10euclidean
  riem10 = Riemann::wrap.euclidean(data_sphere)
  
  ## 12landmark
  data("gorilla", package="Riemann")
  riem12 = Riemann::wrap.landmark(gorilla$male)
  
  ## check altogether
  expect_equal(class(riem01), "riemdata")
  expect_equal(class(riem02), "riemdata")
  expect_equal(class(riem03), "riemdata")
  expect_equal(class(riem04), "riemdata")
  expect_equal(class(riem05), "riemdata")
  expect_equal(class(riem06), "riemdata")
  expect_equal(class(riem07), "riemdata")
  expect_equal(class(riem08), "riemdata")
  expect_equal(class(riem10), "riemdata")
  expect_equal(class(riem12), "riemdata")
})
