#' Data : Left Hands
#' 
#' This dataset contains 10 shapes of 4 subjects's left hands where each shape is represented 
#' by 56 landmark points. For each person, first six shapes are equally spaced sequence 
#' from maximally to minimally spread fingures. The rest are arbitrarily chosen 
#' with two constraints; (1) the palm should face the support and (2) the contour 
#' should contain no crossins.
#' 
#' @usage data(hands)
#' 
#' @examples
#' \donttest{
#' ## LOAD THE DATA 
#' data(hands)
#' 
#' ## VISUALIZE 6 HANDS OF PERSON 1
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' for (i in 1:6){
#'   xx = hands$data[,1,i]
#'   yy = hands$data[,2,i]
#'   plot(xx,yy,"b", cex=0.9)
#' }
#' par(opar)
#' }
#' 
#' @format a named list containing\describe{
#' \item{data}{an \eqn{(56\times 2\times 40)} array of landmarks for 40 subjects.}
#' \item{person}{a length-\eqn{40} vector of subject indices.}
#' }
#' 
#' @references 
#' Stegmann M, Gomez D (2002) "A Brief Introduction to Statistical Shape Analysis." \emph{Informatics and Mathematical Modelling, Technical University of Denmark, DTU.}
#' 
#' @seealso \code{\link{wrap.landmark}}
#' @concept data
"hands"