#' Hierarchical Agglomerative Clustering
#' 
#' Given \eqn{N} observations \eqn{X_1, X_2, \ldots, X_M \in \mathcal{M}}, 
#' perform hierarchical agglomerative clustering with 
#' \pkg{fastcluster} package's implementation.
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param geometry (case-insensitive) name of geometry; either geodesic (\code{"intrinsic"}) or embedded (\code{"extrinsic"}) geometry.
#' @param method agglomeration method to be used. This must be one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
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
#' mydata = list()
#' for (i in 1:10){
#'   tgt = c(1, stats::rnorm(2, sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 11:20){
#'   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' for (i in 21:30){
#'   tgt = c(stats::rnorm(2, sd=0.1), 1)
#'   mydata[[i]] = tgt/sqrt(sum(tgt^2))
#' }
#' myriem = wrap.sphere(mydata)
#' 
#' ## COMPUTE SINGLE AND COMPLETE LINKAGE
#' hc.sing <- riem.hclust(myriem, method="single")
#' hc.comp <- riem.hclust(myriem, method="complete")
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(hc.sing, main="single linkage")
#' plot(hc.comp, main="complete linkage")
#' par(opar)
#' 
#' @seealso  \code{\link[fastcluster]{hclust}}
#' @concept clustering
#' @export
riem.hclust <- function(riemobj, geometry=c("intrinsic","extrinsic"),
                        method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                   "centroid", "median"), members=NULL){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.hclust : input ",DNAME," should be an object of 'riemdata' class."))
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  mymethod = ifelse(missing(method),"single",
                    match.arg(tolower(method), 
                              c("single", "complete", "average", "mcquitty", 
                                "ward.D", "ward.D2", "centroid", "median")))
  mymembers = members
  
  ## COMPUTE DISTANCE, HCLUST, AND RETURN
  pdmat   = stats::as.dist(basic_pdist(riemobj$name, riemobj$data, mygeom))
  fimport = utils::getFromNamespace("hidden_hclust", "maotai")
  hcout   = fimport(pdmat, mymethod, mymembers)
  return(hcout)
}