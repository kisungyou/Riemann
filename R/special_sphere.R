## SOME SPECIAL FUNCTIONS ON SPHERE
#  (01) sphere.runif
#  (02) sphere.utest

# (01) sphere.runif ============================================================
#' Generate Uniform Samples on Sphere
#' 
#' It generates \eqn{n} random samples from \eqn{\mathcal{S}^{p-1}}. For convenient 
#' usage of users, we provide a number of options in terms of the return type.
#' 
#' @param n number of samples to be generated.
#' @param p original dimension (of the ambient space).
#' @param type return type; \describe{
#' \item{\code{"list"}}{a length-\eqn{n} list of length-\eqn{p} vectors.}
#' \item{\code{"matrix"}}{a \eqn{(n\times p)} where rows are unit vectors.}
#' \item{\code{"riemdata"}}{a S3 object. See \code{\link{wrap.sphere}} for more details (\emph{Default}).}
#' }
#' 
#' @return an object from one of the above by \code{type} option.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                       Draw Samples on Sphere
#' #
#' # Multiple return types on S^4 in R^5
#' #-------------------------------------------------------------------
#' dat.list = sphere.runif(n=10, p=5, type="list")
#' dat.matx = sphere.runif(n=10, p=5, type="matrix")
#' dat.riem = sphere.runif(n=10, p=5, type="riemdata")
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept sphere
#' @export
sphere.runif <- function(n, p, type=c("list","matrix","riemdata")){
  # PREPROCESSING
  # parameters
  N = round(n)
  p = round(p) # S^{p-1} in R^p
  
  # return object type
  if (missing(type)){
    retype = "riemdata"
  } else {
    retype = match.arg(tolower(type), c("list","matrix","riemdata"))
  }
  
  # GENERATION
  outmat = runif_sphere(N,p) 
  if (all(retype=="matrix")){
    return(outmat)
  } else {
    outobj = wrap.sphere(outmat)
    if (all(retype=="list")){
      return(outobj$data)
    } else {
      return(outobj)
    }
  }
}

# (02) sphere.utest ============================================================
#' Test of Uniformity on Sphere
#' 
#' Given \eqn{N} observations \eqn{\lbrace X_1, X_2, \ldots, X_M \brace} on 
#' \eqn{\mathcal{S}^{p-1}}, it tests whether the data is distributed uniformly 
#' on the sphere. 
#' 
#' @param spobj a S3 \code{"riemdata"} class for \eqn{N} Sphere-valued data.
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
#' 
#' @examples
#' #-------------------------------------------------------------------
#' #   Compare Rayleigh's original and modified versions of the test
#' #-------------------------------------------------------------------
#' #  Data Generation
#' myobj = sphere.runif(n=100, p=5, type="riemdata")
#' 
#' #  Compare 2 versions : Original vs Modified Rayleigh
#' sphere.utest(myobj, method="rayleigh")
#' sphere.utest(myobj, method="rayleighm")
#' 
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept sphere
#' @export
sphere.utest <- function(spobj, method=c("Rayleigh","RayleighM")){
  ## CHECK INPUT
  check_inputmfd(spobj, "sphere.utest")
  DNAME      = deparse(substitute(spobj)) # borrowed from HDtest
  method.now = ifelse(missing(method),"rayleigh",
                      match.arg(tolower(method),
                                c("rayleigh","rayleighm")))
  
  ## COMPUTE
  # prepare data in 3d array
  data3d <- aux_rvec2array3d(spobj)
  output <- switch(method.now,
                   rayleigh  = sp.utest.Rayleigh(data3d, DNAME, is.original = TRUE),
                   rayleighm = sp.utest.Rayleigh(data3d, DNAME, is.original = FALSE))
  return(output)
}
#' @keywords internal
#' @noRd
sp.utest.Rayleigh <- function(x, dname, is.original=TRUE){
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
  
  if (is.original){
    hname   = "Rayleigh Test of Uniformity on Sphere"
    thestat = S
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  } else {
    hname   = "Modified Rayleigh Test of Uniformity on Sphere"
    term1   = 1/(2*n)
    term2   = 1 - (S/((p*r) + 2))
    thestat = S*(1 - term1*term2)
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  }
  
  # COMPUTATION : DETERMINATION
  Ha      = paste("data is not uniformly distributed on ",p-1,"-sphere.",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}