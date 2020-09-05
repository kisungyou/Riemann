#' Find the Smallest Enclosing Ball
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, find the smallest enclosing ball.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param method (case-insensitive) name of the algorithm to be used as follows;\describe{
#' \item{\code{"aa2013"}}{Arnaudon and Nielsen (2013).}
#' }
#' @param maxiter maximum number of iterations to be run.
#' @param eps tolerance level for stopping criterion.
#' 
#' @return a named list containing \describe{
#' \item{center}{a matrix on \eqn{\mathcal{M}} that minimizes the radius.}
#' \item{radius}{the minimal radius with respect to the \code{center}.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #       Euclidean Example : samples from Standard Normal in R^2
#' #-------------------------------------------------------------------
#' ## GENERATE 25 OBSERVATIONS FROM N(0,I)
#' ndata  = 25
#' mymats = array(0,c(ndata, 2))
#' mydata = list()
#' for (i in 1:ndata){
#'   mydata[[i]] = stats::rnorm(2)
#'   mymats[i,]  = mydata[[i]]
#' }
#' myriem = wrap.euclidean(mydata)
#' 
#' ## COMPUTE
#' sebobj = riem.seb(myriem)
#' center = as.vector(sebobj$center) 
#' radius = sebobj$radius
#' 
#' ## VISUALIZE
#' #  1. prepare the circle for drawing
#' theta  = seq(from=0, to=2*pi, length.out=100)
#' coords = radius*cbind(cos(theta), sin(theta))
#' coords = coords + matrix(rep(center, each=100), ncol=2)
#' 
#' #  2. draw
#' opar <- par(no.readonly=TRUE)
#' par(pty="s")
#' plot(coords, type="l", lwd=2, col="red",
#'      main="Euclidean SEB", xlab="x", ylab="y")
#' points(mymats, pch=19)                           # data
#' points(center[1], center[2], pch=19, col="blue") # center 
#' par(opar)
#' 
#' @references 
#' \insertRef{badoiu_smaller_2003}{Riemann}
#' 
#' \insertRef{arnaudon_approximating_2013}{Riemann}
#' 
#' @concept learning
#' @export
riem.seb <- function(riemobj, method=c("aa2013"), maxiter=50, eps=1e-5){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.seb : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myiter   = max(50, round(maxiter))
  myeps    = as.double(eps)
  mymethod = ifelse(missing(method), "aa2013",
                    match.arg(tolower(method),
                              c("aa2013")))
  
  ## COMPUTE AND RETURN
  output = learning_seb(riemobj$name, riemobj$data, myiter, myeps, mymethod)
  return(output)
}