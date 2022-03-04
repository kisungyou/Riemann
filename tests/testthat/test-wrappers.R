library(Riemann)

test_that("wrappers work", {
  
  ## 01sphere 
  data_sphere = array(0,c(5,3))
  list_sphere = list()
  for (i in 1:5){
    tgt = stats::rnorm(3)
    tgt = tgt/sqrt(sum(tgt^2))
    data_sphere[i,]  = tgt
    list_sphere[[i]] = tgt
  }
  riem01 = Riemann::wrap.sphere(data_sphere)
  list01 = Riemann::wrap.sphere(list_sphere)
  
  ## 02spd
  data_spd = array(0,c(3,3,5))
  list_spd = list()
  for (i in 1:5){
    dat = stats::cov(matrix(rnorm(10*3), ncol=3))
    data_spd[,,i] = dat
    list_spd[[i]] = dat
  }
  riem02 = Riemann::wrap.spd(data_spd)
  list02 = Riemann::wrap.spd(list_spd)
  
  ## 03correlation
  data_corr = array(0,c(3,3,5))
  list_corr = list()
  for (i in 1:5){
    dat = stats::cor(matrix(rnorm(10*3), ncol=3))
    data_corr[,,i] = dat
    list_corr[[i]] = dat 
  }
  riem03 = Riemann::wrap.correlation(data_corr)
  list03 = Riemann::wrap.correlation(list_corr)
  
  ## 04stiefel & 05grassmann
  data_stiefel = array(0,c(4,2,5))
  list_stiefel = list()
  for (i in 1:5){
    dat = qr.Q(qr(matrix(rnorm(4*2), ncol=2)))
    data_stiefel[,,i] = dat 
    list_stiefel[[i]] = dat 
  }
  riem04 = Riemann::wrap.stiefel(data_stiefel)
  list04 = Riemann::wrap.stiefel(list_stiefel)
  riem05 = Riemann::wrap.grassmann(data_stiefel)
  list05 = Riemann::wrap.grassmann(list_stiefel)
  
  ## 06rotation
  data_rotation = array(0,c(3,3,5))
  list_rotation = list()
  for (i in 1:5){
    dat = qr.Q(qr(matrix(rnorm(9), ncol=3)))
    data_rotation[,,i] = dat
    list_rotation[[i]] = dat 
  }
  riem06 = Riemann::wrap.rotation(data_rotation)
  list06 = Riemann::wrap.rotation(list_rotation)
  
  ## 07multinomial
  data_multinomial = array(0,c(5,3))
  list_multinomial = list()
  for (i in 1:5){
    dat = abs(stats::rnorm(3))
    dat = dat/base::sum(dat)
    data_multinomial[i,]  = dat
    list_multinomial[[i]] = dat 
  }
  riem07 = Riemann::wrap.multinomial(data_multinomial)
  list07 = Riemann::wrap.multinomial(list_multinomial)
  
  ## 08spdk
  data_spdk = array(0,c(5,5,3))
  list_spdk = list()
  for (i in 1:3){
    dat = stats::cov(matrix(rnorm(10*5), ncol=5))
    data_spdk[,,i] = dat
    list_spdk[[i]] = dat 
  }
  riem08 = Riemann::wrap.spdk(data_spdk, k=3)
  list08 = Riemann::wrap.spdk(list_spdk, k=3)
  
  ## 10euclidean
  riem10 = Riemann::wrap.euclidean(data_sphere)
  list10 = Riemann::wrap.euclidean(list_sphere)
  
  ## 12landmark
  data("gorilla", package="Riemann")
  list_gorilla = list()
  for (i in 1:dim(gorilla$male)[3]){
    list_gorilla[[i]] = gorilla$male[,,i]
  }
  riem12 = Riemann::wrap.landmark(gorilla$male)
  list12 = Riemann::wrap.landmark(list_gorilla)

  ## check altogether
  expect_equal(class(riem01), "riemdata")
  expect_equal(class(list01), "riemdata")
  
  expect_equal(class(riem02), "riemdata")
  expect_equal(class(list02), "riemdata")
  
  expect_equal(class(riem03), "riemdata")
  expect_equal(class(list03), "riemdata")
  
  expect_equal(class(riem04), "riemdata")
  expect_equal(class(list04), "riemdata")
  
  expect_equal(class(riem05), "riemdata")
  expect_equal(class(list05), "riemdata")
  
  expect_equal(class(riem06), "riemdata")
  expect_equal(class(list06), "riemdata")
  
  expect_equal(class(riem07), "riemdata")
  expect_equal(class(list07), "riemdata")
  
  expect_equal(class(riem08), "riemdata")
  expect_equal(class(list08), "riemdata")
  
  expect_equal(class(riem10), "riemdata")
  expect_equal(class(list10), "riemdata")
  
  expect_equal(class(riem12), "riemdata")
  expect_equal(class(list12), "riemdata")
})
