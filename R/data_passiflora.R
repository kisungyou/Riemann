#' Data : Passiflora Leaves
#' 
#' Passiflora is a genus of about 550 species of flowering plants. This dataset contains 
#' 15 landmarks in 2 dimension of 3319 leaves of 40 species. Papers listed in the 
#' reference section analyzed the data and found 7 clusters.
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
#' plot(pga2d, col=passiflora$class, pch=19, cex=0.7,
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

# library(ggplot2)
# plotter = data.frame(x=pga2d[,1], y=pga2d[,2], cluster=passiflora$class)
# leaves  = ggplot(data=plotter, aes(x=x,y=y,color=cluster)) +
#   geom_point(size=1.5) +
#   theme_bw() +
#   ggtitle("Principal Geodesic Analysis for Passiflora Leaves") +
#   xlab("Dimension 1") +
#   ylab("Dimension 2")
# plot(leaves)
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


# # EXTRINSIC : PLANAR SHAPE IS SAME AS WHAT I HAD
# mydat = passiflora$data[,,1:3]
# myriem = wrap.landmark(mydat)
# riem.pdist(myriem, geometry="intrinsic")
# riem.pdist(myriem, geometry="extrinsic")
# test = array(0,c(3,3))
# for (i in 1:3){
#   datx = myriem$data[[i]]
#   xnow = base::complex(real=datx[,1], imaginary = datx[,2])
#   for (j in 1:3){
#     daty = myriem$data[[j]]
#     ynow = base::complex(real=daty[,1], imaginary = daty[,2])
#     znow = xnow-ynow
#     test[j,i] <- test[i,j] <- sqrt(base::Re(sum(diag(outer(znow,Conj(znow)))))) #== sum(znow*Conj(znow))
#   }
# }
# test
