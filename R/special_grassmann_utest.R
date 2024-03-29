#' Test of Uniformity on Grassmann Manifold
#' 
#' Given the data on Grassmann manifold \eqn{Gr(k,p)}, it tests whether the 
#' data is distributed uniformly.
#' 
#' @param grobj a S3 \code{"riemdata"} class of Grassmann-valued data.
#' @param method (case-insensitive) name of the test method containing \describe{
#' \item{\code{"Bing"}}{Bingham statistic.}
#' \item{\code{"BingM"}}{modified Bingham statistic with better order of error.}
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
#' @seealso \code{\link{wrap.grassmann}}
#' @examples 
#' #-------------------------------------------------------------------
#' #   Compare Bingham's original and modified versions of the test
#' # 
#' # Test 1. sample uniformly from Gr(2,4)
#' # Test 2. use perturbed principal components from 'iris' data in R^4
#' #         which is concentrated around a point to reject H0.
#' #-------------------------------------------------------------------
#' ## Data Generation
#' #  1. uniform data
#' myobj1 = grassmann.runif(n=100, k=2, p=4)
#' 
#' #  2. perturbed principal components
#' data(iris)
#' irdat = list()
#' for (n in 1:100){
#'    tmpdata    = iris[1:50,1:4] + matrix(rnorm(50*4,sd=0.5),ncol=4)
#'    irdat[[n]] = eigen(cov(tmpdata))$vectors[,1:2]
#' }
#' myobj2 = wrap.grassmann(irdat)
#' 
#' ## Test 1 : uniform data
#' grassmann.utest(myobj1, method="Bing")
#' grassmann.utest(myobj1, method="BingM")
#' 
#' ## Tests : iris data
#' grassmann.utest(myobj2, method="bINg")   # method names are 
#' grassmann.utest(myobj2, method="BiNgM")  # CASE - INSENSITIVE !
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{Riemann}
#' 
#' \insertRef{mardia_directional_1999}{Riemann}
#' 
#' @concept grassmann
#' @export
grassmann.utest <- function(grobj, method=c("Bing","BingM")){
  #  CHECK INPUT
  DNAME    = deparse(substitute(grobj)) # borrowed from HDtest
  FNAME    = "grassmann.utest"
  check_inputmfd(grobj, FNAME)
  mymethod = ifelse(missing(method),"bing",
                    match.arg(tolower(method),c("bing","bingm")))
  
  #  COMPUTE AND RETURN
  output <- switch(mymethod,
                   bing  = gr.utest.bing(grobj, DNAME, is.modified = FALSE),
                   bingm = gr.utest.bing(grobj, DNAME, is.modified = TRUE))
  return(output)
}
#' @keywords internal
#' @noRd
gr.utest.bing <- function(grobj, dname, is.modified=TRUE){
  # PREPARE
  p = grobj$size[1] # parameters
  r = grobj$size[2]
  n = length(grobj$data)
  
  # COMPUTE THE MEAN
  Ybar = array(0,c(p,p))
  for (i in 1:n){
    Ytgt = grobj$data[[i]]
    Ybar = Ybar + (Ytgt%*%t(Ytgt))/n
  }
  
  # COMPUTE THE STATISTIC
  S = (((p-1)*p*(p+2))/(2*r*(p-r)))*n*(sum(diag(Ybar%*%Ybar)) - ((r^2)/p))
  
  # BRANCHING
  if (is.modified){
    p2 = p*p
    B0 = ((2*p2*(p-1)*(p+2))-(r*(p-r)*(5*p2 + 2*p + 8)))/(6*r*(p-r)*(p-2)*(p+4))
    B1 = (-(4*p2*(p-1)*(p+2))+(r*(p-r)*((13*p2)+(10*p)-8)))/(3*r*(p-r)*(p2+p+2)*(p-2)*(p+4))
    B2 = (4*((p-(2*r))^2)*(p2+p-2))/(3*r*(p-r)*(p-2)*(p+4)*(p2+p+2)*(p2+p+6))
    
    thestat = S*(1-(1/n)*(B0 + B1*S + B2*(S^2)))
    thedf   = round((p-1)*(p+2)/2)
    hname   = "Modified Bingham Test of Uniformity on Grassmann Manifold"
    pvalue  = stats::pchisq(thestat, df=thedf, lower.tail=FALSE)
  } else {
    hname   = "Bingham Test of Uniformity on Grassmann Manifold"
    thestat = S
    thedf   = round((p-1)*(p+2)/2)
    pvalue  = stats::pchisq(thestat, df=thedf, lower.tail=FALSE)
  }
  
  # DETERMINATION
  Ha      = paste("data is not uniformly distributed on Gr(",r,",",p,").",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = dname)
  class(res) = "htest"
  return(res)
}