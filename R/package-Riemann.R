#' Riemann package
#' 
#' We provide a variety of algorithms for manifold-valued data, including Frechet summaries, 
#' hypothesis testing, clustering, visualization, and other learning tasks. 
#' 
#' @docType package
#' @name Package Introduction
#' @noRd
#' @aliases Riemann-package
#' @import Rdpack
#' @import maotai
#' @import DEoptim
#' @importFrom Matrix nearPD
#' @importFrom CVXR Variable Minimize matrix_trace Problem solve
#' @importFrom lpSolve lp
#' @importFrom T4cluster sc05Z scNJW scSM scUL
#' @importFrom Rdimtools aux.shortestpath
#' @importFrom RiemBase riemfactory
#' @importFrom utils packageVersion getFromNamespace tail
#' @importFrom stats cor rnorm pchisq cov cutree as.dist rnorm runif optimize integrate var kmeans
#' @importFrom Rcpp evalCpp
#' @useDynLib Riemann
NULL
# pack <- "Riemann"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))