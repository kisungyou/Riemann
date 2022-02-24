# methods for S3 class 'mixriem'
#
# () loglkd  : compute the log-likelihood
# () label   : predict the labels 
# () density : evaluate the density



#' Compute the log-likelihood
#' 
#' Given a fitted mixture model \eqn{f(x)} and observations \eqn{x_1, \ldots, x_n \in \calM}, compute the log-likelihood
#' \deqn{L = \log \prod_{i=1}^n f(x_i) = \sum_{i=1}^n \log f(x_i)}.
#' 
#' @param x a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return the log-likelihood.
#' 
#' @export
loglkd <- function(x, newdata){
  UseMethod("loglkd")
}

#' Predict labels of given data for a fitted mixture model
#' 
#' Given a fitted mixture model of \eqn{K} components, predict labels of 
#' observations accordingly.
#' 
#' @param x a fitted mixture model of \code{riemmix} class.
#' @param newdata data of \eqn{n} objects (vectors, matrices) that can be wrapped by one of \code{wrap.*} functions in the \pkg{Riemann} package.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @export
label <- function(x, newdata){
  UseMethod("label")
}

#' @export
density <- function(x, newdata){
  UseMethod("density")
}