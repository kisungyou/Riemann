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
#' @importFrom RiemBase riemfactory
#' @importFrom utils packageVersion getFromNamespace
#' @importFrom stats cor rnorm pchisq cov cutree
#' @importFrom Rcpp evalCpp
#' @useDynLib Riemann
NULL
# pack <- "Riemann"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))