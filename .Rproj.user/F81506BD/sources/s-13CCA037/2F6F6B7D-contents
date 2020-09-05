#' Nonlinear Mean Shift
#' 
#' Given \eqn{N} observations  \eqn{X_1, X_2, \ldots, X_N \in \mathcal{M}}, 
#' perform clustering of the data based on the nonlinear mean shift algorithm. 
#' Gaussian kernel is used with the bandwidth \eqn{h} as of 
#' \deqn{G(x_i, x_j) \propto \exp \left( - \frac{\rho^2 (x_i,x_j)}{h^2} \right)}
#' where \eqn{\rho(x,y)} is geodesic distance between two points \eqn{x,y\in\mathcal{M}}. 
#' Numerically, some of the limiting points that collapse into the same cluster are 
#' not exact. For such purpose, we require \code{maxk} parameter to search the 
#' optimal number of clusters based on \eqn{k}-medoids clustering algorithm 
#' in conjunction with silhouette criterion.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param h bandwidth parameter. The larger the \eqn{h} is, the more blurring is applied.
#' @param maxk maximum number of clusters to determine the optimal number of clusters.
#' @param maxiter maximum number of iterations to be run.
#' @param eps tolerance level for stopping criterion.
#' 
#' @return a named list containing\describe{
#' \item{distance}{an \eqn{(N\times N)} distance between modes corresponding to each data point.}
#' \item{cluster}{a length-\eqn{N} vector of class labels.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #          Example on Sphere : a dataset with three types
#' #
#' # class 1 : 10 perturbed data points near (1,0,0) on S^2 in R^3
#' # class 2 : 10 perturbed data points near (0,1,0) on S^2 in R^3
#' # class 3 : 10 perturbed data points near (0,0,1) on S^2 in R^3
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' set.seed(496)
#' ndata  = 10
#' mydata = list()
#' for (i in 1:ndata){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in (ndata+1):(2*ndata)){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in ((2*ndata)+1):(3*ndata)){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' mylabs = rep(c(1,2,3), each=ndata)
#' 
#' ## RUN NONLINEAR MEANSHIFT FOR DIFFERENT 'h' VALUES
#' run1 = riem.nmshift(myriem, maxk=10, h=0.1)
#' run2 = riem.nmshift(myriem, maxk=10, h=1)
#' run3 = riem.nmshift(myriem, maxk=10, h=10)
#' 
#' ## MDS FOR VISUALIZATION
#' mds2d = riem.mds(myriem, ndim=2)$embed
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' plot(mds2d, pch=19, main="label : h=0.1", col=run1$cluster)
#' plot(mds2d, pch=19, main="label : h=1",   col=run2$cluster)
#' plot(mds2d, pch=19, main="label : h=10",  col=run3$cluster)
#' image(run1$distance[,30:1], axes=FALSE, main="distance : h=0.1")
#' image(run2$distance[,30:1], axes=FALSE, main="distance : h=1")
#' image(run3$distance[,30:1], axes=FALSE, main="distance : h=10")
#' par(opar)
#' 
#' @references 
#' \insertRef{subbarao_nonlinear_2009}{Riemann}
#' 
#' @concept clustering
#' @export
riem.nmshift <- function(riemobj, h=1, maxk=5, maxiter=50, eps=1e-5){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.nmshift : input ",DNAME," should be an object of 'riemdata' class."))
  }
  myh    = max(0, as.double(h))    # min(max(as.double(h),0),1e-5)
  myiter = max(50, round(maxiter))
  myeps  = min(max(as.double(eps),0),1e-5)
  
  ## COMPUTE LIMIT POINTS AND PAIRWISE DISTANCE
  tmprun  = clustering_nmshift(riemobj$name, riemobj$data, myh, myiter, myeps)
  limit3d = tmprun$limits
  pdmat   = tmprun$distance
  objmat  = stats::as.dist(pdmat)
  
  ## USE K-MEDOIDS + SILHOUETTE FOR DETERMINATION
  mink = 2
  maxk = max(min(round(maxk), nrow(pdmat)-1), (mink+1))
  
  fimport = utils::getFromNamespace("hidden_kmedoids_best", "maotai")
  hcout   = fimport(objmat, mink, maxk)
  k.opt   = as.integer(hcout$opt.k)
  k.label = hcout$label[,(k.opt-1)]
  
  # ## USE HCLUST + SINGLE LINKAGE + SILHOUETTE FOR DETERMINATION
  # #  prepare for multiple k inspection
  # mink = 2
  # maxk = max(min(round(maxk), nrow(pdmat)-1), (mink+1))
  # veck = seq(from=mink, to=maxk, by=1)
  # silk = rep(0, length(veck))
  # 
  # #  prepare for hclust
  # fimport = utils::getFromNamespace("hidden_hclust", "maotai")
  # hcout   = fimport(objmat, "single", NULL)
  # 
  # #  compute for multiple k's
  # silfunc = utils::getFromNamespace("hidden_silhouette", "maotai")
  # for (i in 1:length(veck)){
  #   labels = stats::cutree(hcout, k=round(veck[i]))
  #   silrun = silfunc(objmat, labels)
  #   silk[i] = silrun$global
  # }
  # k.opt    = veck[which.max(silk)]
  # k.label  = as.integer(as.factor(stats::cutree(hcout, k=k.opt)))
  # dat.nrow = dim(limit3d)[1]
  # dat.ncol = dim(limit3d)[2]

  # ## COMPUTE MEAN PER CLASS
  # out.means = array(0,c(dat.nrow, dat.ncol, k.opt))
  # for (i in 1:k.opt){
  #   idk = which(k.label==i)
  #   if (length(idk)==1){
  #     out.means[,,i] = limit3d[,,idk]
  #   } else {
  #     # heuristic : choose the median
  #     partdist = pdmat[idk,idk]
  #     idmin    = idk[which.min(base::rowSums(partdist^2))]
  #     out.means[,,i] = (limit3d[,,idmin])
  # 
  #     # tmplist = list()
  #     # for (j in 1:length(idk)){
  #     #   tmplist[[j]] = wrap_vec2mat(limit3d[,,idk[j]])
  #     # }
  #     # tmpweight = rep(1/length(idk), length(idk))
  #     # tmpmean   = inference_mean_intrinsic(riemobj$name, tmplist, tmpweight, myiter, myeps)
  #     # out.means[,,i] = wrap_vec2mat(tmpmean$mean)
  #   }
  # }
  
  ## WRAP AND RETURN
  output = list()
  output$distance = pdmat
  output$cluster  = k.label
  return(output)
}

# p=2
# x = c(1,rep(0,p))
# y = -x
# acos(sum(x*y))

# which(run3$clustering==9)
# partdata = mydata[c(24,30)]
# partriem = wrap.sphere(partdata)
# riem.mean(partriem)
