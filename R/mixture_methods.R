# methods for S3 class 'mixriem'
# check with --as-cran --run-donttest --use-valgrind
#
# () loglkd  : compute the log-likelihood
# () label   : predict the labels 
# () density : evaluate the density



#' Compute the log-likelihood
#' 
#' Given a fitted mixture model \eqn{f(x)} and observations \eqn{x_1, \ldots, x_n \in \calM}, compute the log-likelihood
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
#' #                 FITTING THE MODEL
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit the model with different numbers of clusters
#' k2 = moSN(locations, k=2)
#' k3 = moSN(locations, k=3)
#' k4 = moSN(locations, k=4)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed2, col=k2$cluster, pch=19, main="K=2")
#' plot(embed2, col=k3$cluster, pch=19, main="K=3")
#' plot(embed2, col=k4$cluster, pch=19, main="K=4")
#' par(opar)
#' 
#' # ---------------------------------------------------- #
#' #                   USE S3 METHODS
#' # ---------------------------------------------------- #
#' # Use the same 'locations' data as new data 
#' # (1) log-likelihood
#' newloglkd = round(loglkd(k3, locations), 3)
#' print(paste0("Log-likelihood for K=3 model fit : ", newloglkd))
#' 
#' # (2) label
#' newlabel = label(k3, locations)
#' 
#' # (3) density
#' newdensity = density(k3, locations)
#' }
#' 
#' @export
loglkd <- function(object, newdata){
  UseMethod("loglkd")
}

#' Predict labels of given data for a fitted mixture model
#' 
#' Given a fitted mixture model of \eqn{K} components, predict labels of 
#' observations accordingly.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------- #
#' #                 FITTING THE MODEL
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit the model with different numbers of clusters
#' k2 = moSN(locations, k=2)
#' k3 = moSN(locations, k=3)
#' k4 = moSN(locations, k=4)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed2, col=k2$cluster, pch=19, main="K=2")
#' plot(embed2, col=k3$cluster, pch=19, main="K=3")
#' plot(embed2, col=k4$cluster, pch=19, main="K=4")
#' par(opar)
#' 
#' # ---------------------------------------------------- #
#' #                   USE S3 METHODS
#' # ---------------------------------------------------- #
#' # Use the same 'locations' data as new data 
#' # (1) log-likelihood
#' newloglkd = round(loglkd(k3, locations), 3)
#' print(paste0("Log-likelihood for K=3 model fit : ", newloglkd))
#' 
#' # (2) label
#' newlabel = label(k3, locations)
#' 
#' # (3) density
#' newdensity = density(k3, locations)
#' }
#' 
#' @param object a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @export
label <- function(object, newdata){
  UseMethod("label")
}



#' Evaluate density of given data for a fitted model
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
#' #                 FITTING THE MODEL
#' # ---------------------------------------------------- #
#' # Load the 'city' data and wrap as 'riemobj'
#' data(cities)
#' locations = cities$cartesian
#' embed2    = array(0,c(60,2)) 
#' for (i in 1:60){
#'    embed2[i,] = sphere.xyz2geo(locations[i,])
#' }
#' 
#' # Fit the model with different numbers of clusters
#' k2 = moSN(locations, k=2)
#' k3 = moSN(locations, k=3)
#' k4 = moSN(locations, k=4)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed2, col=k2$cluster, pch=19, main="K=2")
#' plot(embed2, col=k3$cluster, pch=19, main="K=3")
#' plot(embed2, col=k4$cluster, pch=19, main="K=4")
#' par(opar)
#' 
#' # ---------------------------------------------------- #
#' #                   USE S3 METHODS
#' # ---------------------------------------------------- #
#' # Use the same 'locations' data as new data 
#' # (1) log-likelihood
#' newloglkd = round(loglkd(k3, locations), 3)
#' print(paste0("Log-likelihood for K=3 model fit : ", newloglkd))
#' 
#' # (2) label
#' newlabel = label(k3, locations)
#' 
#' # (3) density
#' newdensity = density(k3, locations)
#' }
#' 
#' @export
density <- function(object, newdata){
  UseMethod("density")
}