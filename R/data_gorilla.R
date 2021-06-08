#' Data : Gorilla Skull
#' 
#' This is 29 male and 30 female gorillas' skull landmark data where each 
#' individual is represented as 8-ad/landmarks in 2 dimensions. This is a 
#' re-arranged version of the data from \pkg{shapes} package.
#' 
#' @usage data(gorilla)
#' 
#' @format a named list containing\describe{
#' \item{male}{a 3d array of size \eqn{(8\times 2\times 29)}}
#' \item{female}{a 3d array of size \eqn{(8\times 2\times 30)}}
#' }
#' 
#' @examples 
#' data(gorilla)                               # load the data
#' riem.female = wrap.landmark(gorilla$female) # wrap as RIEMOBJ
#' opar <- par(no.readonly=TRUE)
#' for (i in 1:30){
#'   if (i < 2){
#'     plot(riem.female$data[[i]], cex=0.5, 
#'          xlim=c(-1,1)/2, ylim=c(-2,2)/5,
#'          main="30 female gorilla skull preshapes",
#'          xlab="dimension 1", ylab="dimension 2")
#'     lines(riem.female$data[[i]])
#'   } else {
#'     points(riem.female$data[[i]], cex=0.5)
#'     lines(riem.female$data[[i]])
#'   }
#' }
#' par(opar)
#' 
#' @references 
#' \insertRef{dryden_statistical_2016}{Riemann}
#' 
#' Reno PL, Meindl RS, McCollum MA, Lovejoy CO (2003). "Sexual dimorphism in Australopithecus afarensis was similar to that of modern humans." \emph{Proceedings of the National Academy of Sciences}, 100(16):9404–9409.
#' 
#' @seealso \code{\link{wrap.landmark}}
#' @concept data
"gorilla"