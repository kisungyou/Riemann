# methods for S3 class 'mixriem'
# check with --as-cran --run-donttest --use-valgrind
#
# () loglkd  : compute the log-likelihood
# () label   : predict the labels 
# () density : evaluate the density



#' S3 method for mixture model : log-likelihood
#' 
#' Given a fitted mixture model \eqn{f(x)} and observations \eqn{x_1, \ldots, x_n \in \mathcal{M}}, compute the log-likelihood
#' \deqn{L = \log \prod_{i=1}^n f(x_i) = \sum_{i=1}^n \log f(x_i)}.
#' 
#' @param object a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return the log-likelihood.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------- #
#' #            FIT A MODEL & APPLY THE METHOD
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit a model
#' k3 = moSN(locations, k=3)
#' 
#' # Evaluate
#' newloglkd = round(loglkd(k3, locations), 3)
#' print(paste0("Log-likelihood for K=3 model fit : ", newloglkd))
#' }
#' 
#' @concept utility
#' @export
loglkd <- function(object, newdata){
  UseMethod("loglkd")
}

#' S3 method for mixture model : predict labels
#' 
#' Given a fitted mixture model of \eqn{K} components, predict labels of 
#' observations accordingly.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------- #
#' #            FIT A MODEL & APPLY THE METHOD
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit a model
#' k3 = moSN(locations, k=3)
#' 
#' # Evaluate
#' newlabel = label(k3, locations)
#' }
#' 
#' @param object a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @concept utility
#' @export
label <- function(object, newdata){
  UseMethod("label")
}



#' S3 method for mixture model : evaluate density
#' 
#' Compute density for a fitted mixture model.
#' 
#' @param object a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------- #
#' #            FIT A MODEL & APPLY THE METHOD
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit a model
#' k3 = moSN(locations, k=3)
#' 
#' # Evaluate 
#' newdensity = density(k3, locations)
#' }
#' 
#' @concept utility
#' @export
density <- function(object, newdata){
  UseMethod("density")
}