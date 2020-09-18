#' Data : Passiflora Leaves
#' 
#' Passiflora is a genus of about 550 species of flowering plants. This dataset contains 
#' 15 landmarks in 2 dimension of 3319 leaves of 40 species. Papers listed in the 
#' reference section analyzed the data and found 7 clusters.
#' 
#' 
#' @usage data(passiflora)
#' 
#' @format a named list containing\describe{
#' \item{data}{a 3d array of size \eqn{(15\times 2\times 3319)}.}
#' \item{species}{a length-\eqn{3319} vector of 40 species factors.}
#' \item{class}{a length-\eqn{3319} vector of 7 cluster factors.}
#' }
#' 
#' @examples 
#' data(passiflora)                         # load the data
#' riemobj = wrap.landmark(passiflora$data) # wrap as RIEMOBJ
#' pga2d   = riem.pga(riemobj)$embed        # embedding via PGA
#' 
#' opar <- par(no.readonly=TRUE)            # visualize
#' plot(pga2d, col=passiflora$class, pch=19, 
#'      main="PGA Embedding of Passiflora Leaves",
#'      xlab="dimension 1", ylab="dimension 2")
#' par(opar)
#' 
#' @references 
#' \insertRef{chitwood_divergent_2017}{Riemann}
#' 
#' \insertRef{chitwood_morphometric_2017}{Riemann}
#' 
#' @seealso \code{\link{wrap.landmark}}
#' @concept data
"passiflora"


# data("passiflora")
# riemobj = wrap.landmark(passiflora$data)
# d2.pga  = riem.pga(riemobj, ndim=2)
# d2.tsne = riem.tsne(riemobj, ndim=2)
# d2.mds  = riem.mds(riemobj, ndim=2)
# 
# dmat = riem.rmml(riemobj, passiflora$class)
# d2.rmml = cmdscale(as.dist(dmat), k = 2)
# 
# par(mfrow=c(1,4))
# plot(d2.pga$embed,  xlab="x", ylab="y", col=passiflora$class, pch=19, main="PGA")
# plot(d2.mds$embed,  xlab="x", ylab="y", col=passiflora$class, pch=19, main="MDS")
# plot(d2.tsne$embed, xlab="x", ylab="y", col=passiflora$class, pch=19, main="t-SNE")
# plot(d2.rmml,       xlab="x", ylab="y", col=passiflora$class, pch=19, main="RMML")
