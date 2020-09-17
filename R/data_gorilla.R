#' Gorilla skull data
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
#' # LOAD AND WRAP AS RIEMOBJ
#' data(gorilla)
#' riem.female = wrap.landmark(gorilla$female)
#' riem.male   = wrap.landmark(gorilla$male)
#' 
#' @references 
#' \insertRef{dryden_statistical_2016}{Riemann}
#' 
#' \insertRef{reno_sexual_2003}{Riemann}
#' 
#' @seealso \code{\link{wrap.landmark}}
#' @concept data
"gorilla"