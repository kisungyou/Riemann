#' Two-Sample Test modified from Biswas and Ghosh (2014)
#' 
#' Given \eqn{M} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}} and 
#' \eqn{N} observations \eqn{Y_1, Y_2, \ldots, Y_N \in \mathcal{M}}, perform the permutation test of equal distribution
#' \deqn{H_0~:~\mathcal{P}_X = \mathcal{P}_Y}
#' by the method from Biswas and Ghosh (2014). The method, originally proposed 
#' for Euclidean-valued data, is adapted to the general Riemannian manifold 
#' with intrinsic/extrinsic distance. 
#' 
#' @param riemobj1 a S3 \code{"riemdata"} class for \eqn{M} manifold-valued data.
#' @param riemobj2 a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param ... extra parameters including\describe{
#' \item{nperm}{the number of permutations (default: 999).}
#' }
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with two types
#' #
#' # class 1 : 20 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 30 perturbed data points near (0,1,0) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' mydata1 = list()
#' mydata2 = list()
#' for (i in 1:20){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata1[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 1:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata2[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem1 = wrap.sphere(mydata1)
#' myriem2 = wrap.sphere(mydata2)
#' 
#' ## PERFORM PERMUTATION TEST
#' #  it is expected to return a very small number.
#' \donttest{
#' riem.test2bg14(myriem1, myriem2, nperm=999)
#' }
#' 
#' \dontrun{
#' ## CHECK WITH EMPIRICAL TYPE-1 ERROR
#' set.seed(777)
#' ntest = 1000
#' pvals = rep(0,ntest)
#' 
#' for (i in 1:ntest){
#'   X = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Y = cbind(matrix(rnorm(30*2, sd=0.1),ncol=2), rep(1,30))
#'   Xnorm = X/sqrt(rowSums(X^2))
#'   Ynorm = Y/sqrt(rowSums(Y^2))
#'   
#'   Xriem = wrap.sphere(Xnorm)
#'   Yriem = wrap.sphere(Ynorm)
#'   pvals[i] = riem.test2bg14(Xriem, Yriem, nperm=999)$p.value
#' }
#' 
#' emperr = round(sum((pvals <= 0.05))/ntest, 5)
#' print(paste0("* EMPIRICAL TYPE-1 ERROR=", emperr))
#' }
#' 
#' @references
#' \insertRef{biswas_nonparametric_2014a}{Riemann}
#' 
#' \insertRef{you_revisiting_2020a}{Riemann}
#' 
#' @concept inference
#' @export
riem.test2bg14 <- function(riemobj1, riemobj2, geometry=c("intrinsic","extrinsic"), ...){
  ## INPUTS : EXPLICIT
  DNAME1 = paste0("'",deparse(substitute(riemobj1)),"'")
  DNAME2 = paste0("'",deparse(substitute(riemobj2)),"'")
  if (!inherits(riemobj1,"riemdata")){
    stop(paste0("* riem.test2bg14 : input ",DNAME1," should be an object of 'riemdata' class."))
  }
  if (!inherits(riemobj2,"riemdata")){
    stop(paste0("* riem.test2bg14 : input ",DNAME2," should be an object of 'riemdata' class."))
  }
  mygeometry = ifelse(missing(geometry),"intrinsic",
                      match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  m = length(riemobj1$data)
  n = length(riemobj2$data)
  
  ## INPUTS : IMPLICIT
  param    = list(...)
  pnames   = names(param)
  mynperm  = ifelse(("nperm"%in%pnames), max(9, round(param$nperm)), 999)
  
  ## COMPUTE PAIRWISE DISTANCE AND PREPARE
  DXY = basic_pdist(riemobj1$name, c(riemobj1$data, riemobj2$data), mygeometry)
  DX0 = DXY[1:m,1:m]                   # under null
  DY0 = DXY[(m+1):(m+n),(m+1):(m+n)]
  DZ0 = DXY[1:m,(m+1):(m+n)]
  Tmn = R_eqdist_2014BG_statistic(DX0,DY0,DZ0)
  
  ## PERMUTATION
  Tvec = rep(0,mynperm)
  for (i in 1:mynperm){
    idx = sample(1:(m+n), m, replace=FALSE)
    idy = setdiff(1:(m+n), idx)
    
    DX1 = DXY[idx,idx]
    DY1 = DXY[idy,idy]
    DZ1 = DXY[idx,idy]
    Tvec[i] = R_eqdist_2014BG_statistic(DX1,DY1,DZ1)
  }
  pvalue = (sum(Tvec>=Tmn)+1)/(mynperm+1)
  
  ## WRAP
  dataname = paste0(DNAME1," and ",DNAME2)
  mfdname  = wrap_mfd2full(riemobj1$name)
  hname    = paste0("Two-Sample Test on ",mfdname," as of Biswas and Ghosh (2014)")
  Ha       = "two distributions are not equal."
  names(Tmn) = "Tmn"
  
  res   = list(statistic=Tmn, p.value=pvalue, alternative = Ha, method=hname, data.name=dataname)
  class(res) = "htest"
  return(res)
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
R_eqdist_2014BG_statistic <- function(DX,DY,DXY){
  m = nrow(DXY)
  n = ncol(DXY)
  
  muff = sum(DX[upper.tri(DX)])/(m*(m-1)/2)
  mufg = sum(DXY)/(m*n)
  mugg = sum(DY[upper.tri(DY)])/(n*(n-1)/2)
  
  vec1 = c(muff,mufg)
  vec2 = c(mufg,mugg)
  output = sum((vec1-vec2)^2)
  return(output)
}