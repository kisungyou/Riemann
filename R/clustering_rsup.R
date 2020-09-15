#' Riemannian Self-Updating Process
#' 
#' @param riemobj a S3 \code{"riemdata"} class for \eqn{N} manifold-valued data.
#' @param weight d
#' @param c d
#' @param maxiter maximum number of iterations to be run.
#' @param eps tolerance level for stopping criterion.
#' 

#' 
#' @return a named list containing\describe{
#' \item{distance}{an \eqn{(N\times N)} distance between modes corresponding to each data point.}
#' \item{cluster}{a length-\eqn{N} vector of class labels.}
#' }
#' @keywords internal
#' @noRd
riem.rsup <- function(riemobj, weight=NULL, c=5, maxiter=50, eps=1e-5){
  ## PREPARE
  N = length(riemobj$data)
  myc    = max(sqrt(.Machine$double.eps), as.double(c))
  mymult = 1/myc
  myiter = max(50, as.integer(maxiter))
  myeps  = min(max(as.double(eps),0),1e-5)
  
  if ((is.null(weight)&&(length(weight)==0))){
    myweight = rep(1/N, N)
  } else {
    myweight = as.double(weight)
    if (length(weight)!=N){
      stop("* riem.rsup : 'weight' vector should be of matching length as the given data.")
    }
    if (any(weight)<=0){
      stop("* riem.rsup : 'weight' vector should contain nonnegative real numbers only.")
    }
  }
  
  ## RUN THE ALGORITHM IN CPP
  tmprun = clustering_sup_intrinsic(riemobj$name, riemobj$data, myweight, mymult, myiter, myeps)
  limit3d = tmprun$limits
  pdmat   = tmprun$distance
  objmat  = stats::as.dist(pdmat)
  
  ## USE K-MEDOIDS + SILHOUETTE FOR DETERMINATION
  mink = 2
  maxk = min(10, round(N/2))
  
  fimport = utils::getFromNamespace("hidden_kmedoids_best", "maotai")
  hcout   = fimport(objmat, mink, maxk)
  k.opt   = as.integer(hcout$opt.k)
  k.label = hcout$label[,(k.opt-1)]
  
  ## WRAP AND RETURN
  output = list()
  output$distance = pdmat
  output$cluster  = k.label
  return(output)
}