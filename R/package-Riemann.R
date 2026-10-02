#' Riemann package
#'
#' We provide a variety of algorithms for manifold-valued data, including Frechet summaries,
#' hypothesis testing, clustering, visualization, and other learning tasks.
#'
#' @name Riemann-package
#' @aliases Riemann-package
#' @section License:
#' Author-owned software in the September 21, 2026 revision of version 0.2.0 is
#' licensed under the GNU General Public License, version 3 (GPL-3).
#' Earlier MIT releases and preserved archives retain their original licenses.
#' Third-party notices are in the installed \file{COPYRIGHTS} file.
#' Dataset attribution and any data-specific terms are separate; consult the
#' dataset help pages. The separate manuscript and replication workspace has
#' its own license notices.
#' @import Rdpack
#' @import maotai
#' @import DEoptim
#' @importFrom Matrix nearPD
#' @importFrom T4transport wassersteinD
#' @importFrom T4cluster sc05Z scNJW scSM scUL
#' @importFrom utils packageVersion getFromNamespace tail
#' @importFrom stats cor rnorm pchisq cov cutree as.dist rnorm runif optimize integrate var kmeans pnorm hclust
#' @importFrom Rcpp evalCpp
#' @useDynLib Riemann
"_PACKAGE"
# pack <- "Riemann"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
