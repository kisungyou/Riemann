## SOME SPECIAL FUNCTIONS ON STIEFEL 
#  (01) stiefel.runif
#  (02) stiefel.utest

#  (01) stiefel.runif ==========================================================
#' Generate Uniform Samples on Stiefel Manifold
#' 
#' It generates \eqn{n} random samples from Stiefel manifold \eqn{St(k,p)}.
#' 
#' @param n number of samples to be generated.
#' @param k dimension of the frame.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of \eqn{(p\times k)} basis of \eqn{k}-frames.}
#' \item{\code{"array"}}{a \eqn{(p\times k\times n)} 3D array whose slices are \eqn{k}-frame basis.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.stiefel}} for more details.}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' @seealso \code{\link{wrap.stiefel}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Draw Samples on Stiefel Manifold 
#' #
#' # Try Different Return Types with 3 Observations of 5-frames in R^10
#' #-------------------------------------------------------------------
#' #  GENERATION
#' dat.list = stiefel.runif(n=3, k=5, p=10, type="list")
#' dat.arr3 = stiefel.runif(n=3, k=5, p=10, type="array")
#' dat.riem = stiefel.runif(n=3, k=5, p=10, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @concept stiefel
#' @export
stiefel.runif <- function(n, k, p, type=c("list","array","riemdata")){
  ## PREPROCESSING
  N = round(n)
  p = round(p) # k-frame in R^p
  k = round(k)
  if (k > p){
    stop("* stiefel.runif : 'k <= p' is a required condition.")
  }
  ## MAIN COMPUTATION
  #  return object type
  retype = ifelse(missing(type), "riemdata",
                  match.arg(tolower(type), c("list","array","riemdata")))
  tmpout = runif_stiefel(p,k,N) # C++ version

  ## RETURN
  if (all(retype=="array")){
    return(tmpout)
  } else {
    tmpobj = wrap.stiefel(tmpout)
    if (all(retype=="list")){
      return(tmpobj$data)
    } else {
      return(tmpobj)
    }
  }
}

#  (02) stiefel.utest ==========================================================
#' Test of Uniformity on Stiefel Manifold
#' 
#' Given the data on Stiefel manifold \eqn{St(k,p)}, it tests whether the 
#' data is distributed uniformly.
#' 
#' @param stobj a S3 \code{"riemdata"} class for \eqn{N} Stiefel-valued data.
#' @param method (case-insensitive) name of the test method containing \describe{
#' \item{\code{"Rayleigh"}}{original Rayleigh statistic.}
#' \item{\code{"RayleighM"}}{modified Rayleigh statistic with better order of error.}
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
#' @seealso \code{\link{wrap.stiefel}}
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #   Compare Rayleigh's original and modified versions of the test
#' # 
#' # Test 1. sample uniformly from St(2,4)
#' # Test 2. use perturbed principal components from 'iris' data in R^4
#' #         which is concentrated around a point to reject H0.
#' #-------------------------------------------------------------------
#' ## DATA GENERATION
#' #  1. uniform data
#' myobj1 = stiefel.runif(n=100, k=2, p=4)
#' 
#' #  2. perturbed principal components
#' data(iris)
#' irdat = list()
#' for (n in 1:100){
#'    tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
#'    irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
#' }
#' myobj2 = wrap.stiefel(irdat)
#' 
#' ## TEST
#' #  1. uniform data
#' stiefel.utest(myobj1, method="Rayleigh")
#' stiefel.utest(myobj1, method="RayleighM")
#' 
#' #  2. concentrated data
#' stiefel.utest(myobj2, method="rayleIgh")   # method names are 
#' stiefel.utest(myobj2, method="raYleiGhM")  # CASE - INSENSITIVE !
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @concept stiefel
#' @export
stiefel.utest <- function(stobj, method=c("Rayleigh","RayleighM")){
  ## CHECK INPUT
  FNAME      = "stiefel.utest"
  DNAME      = deparse(substitute(stobj)) # borrowed from HDtest
  check_inputmfd(stobj, FNAME)
  method.now = ifelse(missing(method),"rayleigh",
                      match.arg(tolower(method),c("rayleigh","rayleighm")))
  
  ## COMPUTATION
  #  prepare data in 3d array
  data3d <- aux_rmat2array3d(stobj)
  output <- switch(method.now,
                   rayleigh  = st.utest.Rayleigh(data3d, DNAME, is.modified = FALSE),
                   rayleighm = st.utest.Rayleigh(data3d, DNAME, is.modified = TRUE))
  return(output)
}
#' @keywords internal
#' @noRd
st.utest.Rayleigh <- function(x, dname, is.modified=FALSE){
  # Take the version from RiemStiefel {will be deprecated}
  p = dim(x)[1]      # r-frame in R^p with n observations
  r = dim(x)[2]
  n = dim(x)[3]
  xbar = array(0,c(p,r))
  for (i in 1:p){
    for (j in 1:r){
      xbar[i,j] = base::mean(as.vector(x[i,j,]))
    }
  }
  S    = p*n*sum(diag(t(xbar)%*%xbar))
  
  if (is.modified){
    hname   = "Modified Rayleigh Test of Uniformity on Stiefel Manifold"
    term1   = 1/(2*n)
    term2   = 1 - (S/((p*r) + 2))
    thestat = S*(1 - term1*term2)
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)  
  } else {
    hname   = "Rayleigh Test of Uniformity on Stiefel Manifold"
    thestat = S
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  }
  
  # COMPUTATION : DETERMINATION
  Ha      = paste("data is not uniformly distributed on St(",r,",",p,").",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}