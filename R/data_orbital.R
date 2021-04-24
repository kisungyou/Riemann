#' Data : Normal Vectors to the Orbital Planes of the 9 Planets
#' 
#' The 9 planets in our solar system are evolving the sun via their own orbits. 
#' This data provides normal vector of the orbital planes. Normal vectors are 
#' unit-norm vectors, so that they are thought to reside on 2-dimensional sphere. 
#' 
#' @usage data(orbital)
#' 
#' @examples
#' ## LOAD THE DATA AND WRAP AS RIEMOBJ
#' data(orbital)
#' myorb = wrap.sphere(orbital)
#' 
#' ## VISUALIZE
#' mds2d = riem.mds(myorb)$embed
#' opar <- par(no.readonly=TRUE)
#' plot(mds2d, main="9 Planets", pch=19, xlab="x", ylab="y")
#' par(opar)
#' 
#' @format an \eqn{(9\times 3)} matrix where each row is a normal vector for a planet.
#' 
#' @seealso \code{\link{wrap.sphere}}
#' @concept data
"orbital"