#' Fréchet Analysis of Variance
#' 
#' Given sets of manifold-valued data \eqn{X^{(1)}_{1:{n_1}}, X^{(2)}_{1:{n_2}}, \ldots, X^{(m)}_{1:{n_m}}}, 
#' performs analysis of variance to test equality of distributions. This means, small \eqn{p}-value implies that 
#' at least one of the equalities does not hold. 
#' 
#' @param ... S3 objects of \code{riemdata} class for manifold-valued data.
#' @param nperm the number of permutations for resampling-based test.
#' @param maxiter maximum number of iterations to be run.
#' @param eps tolerance level for stopping criterion.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #            Example on Sphere : Uniform Samples
#' #
#' #  Each of 4 classes consists of 20 uniform samples from uniform 
#' #  density on 2-dimensional sphere S^2 in R^3.
#' #-------------------------------------------------------------------
#' ## PREPARE DATA OF 4 CLASSES
#' ndata  = 20
#' class1 = list()
#' class2 = list()
#' class3 = list()
#' class4 = list()
#' for (i in 1:ndata){
#'   tmp = matrix(rnorm(4*3), ncol=3)
#'   tmp = tmp/sqrt(rowSums(tmp^2))
#'   
#'   class1[[i]] = tmp[1,]
#'   class2[[i]] = tmp[2,]
#'   class3[[i]] = tmp[3,]
#'   class4[[i]] = tmp[4,]
#' }
#' obj1 = wrap.sphere(class1)
#' obj2 = wrap.sphere(class2)
#' obj3 = wrap.sphere(class3)
#' obj4 = wrap.sphere(class4)
#' 
#' ## RUN THE ASYMPTOTIC TEST
#' riem.fanova(obj1, obj2, obj3, obj4)
#' 
#' \donttest{
#' ## RUN THE PERMUTATION TEST WITH MANY PERMUTATIONS
#' riem.fanovaP(obj1, obj2, obj3, obj4, nperm=9999)
#' }
#' 
#' @references 
#' \insertRef{dubey_frechet_2019}{Riemann}
#' 
#' @name riem.fanova
#' @concept inference
#' @rdname riem.fanova
NULL

#' @rdname riem.fanova
#' @export
riem.fanova <- function(..., maxiter=50, eps=1e-5){
  ## PREPARE
  #  data
  datalist = base::list(...)
  k = length(datalist)
  for (i in 1:k){
    riemobj = datalist[[i]]
    if (!inherits(riemobj, "riemdata")){
      stop(paste0("* riem.fanova : ",i,"-th input should be an object of 'riemdata' class."))
    }
  }
  # hypothesis testing argument
  DNAME = ""
  DOBJ  = as.list(substitute(list(...)))[-1L]
  for (i in 1:(k-1)){
    DNAME = paste0(DNAME, as.character(DOBJ[[i]]), ", ")
  }
  DNAME = paste0(DNAME, "and ",as.character(DOBJ[[k]]))
  MNAME = riemobj$name
  # iteration
  myiter = max(50, round(maxiter))
  myeps  = min(max(as.double(eps),0),1e-5)
  
  ## COMPUTATION
  # Step 1. compute pooled frechet objective
  pooled.data = list()
  for (i in 1:k){
    pooled.data = c(pooled.data, datalist[[i]]$data)
  }
  pooled.ndata   = length(pooled.data)
  pooled.weight  = rep(1/pooled.ndata, pooled.ndata)
  frechet.pooled = inference_mean_intrinsic(MNAME, pooled.data, pooled.weight, myiter, myeps)
  
  # Step 2. compute individual frechet objective
  frechet.each = list()
  for (i in 1:k){
    tgt.ndata  = length(datalist[[i]]$data)
    tgt.weight = rep(1/tgt.ndata, tgt.ndata)
    frechet.each[[i]] = inference_mean_intrinsic(MNAME, datalist[[i]]$data, tgt.weight, myiter, myeps)
  }
  
  # Step 3. get distance information only and compute via 'common_fanova'
  dist.pooled = as.vector(frechet.pooled$distvec)
  dist.class  = list()
  for (i in 1:k){
    dist.class[[i]] = as.vector(frechet.each[[i]]$distvec)
  }
  statinfo = common_fanova(dist.pooled, dist.class, MNAME, DNAME)
  
  # RETURN
  return(statinfo)
}

#' @rdname riem.fanova
#' @export
riem.fanovaP <- function(..., maxiter=50, eps=1e-5, nperm=99){
  ## PREPARE
  #  data
  datalist = base::list(...)
  k = length(datalist)
  for (i in 1:k){
    riemobj = datalist[[i]]
    if (!inherits(riemobj, "riemdata")){
      stop(paste0("* riem.fanovaP : ",i,"-th input should be an object of 'riemdata' class."))
    }
  }
  # hypothesis testing argument
  DNAME = ""
  DOBJ  = as.list(substitute(list(...)))[-1L]
  for (i in 1:(k-1)){
    DNAME = paste0(DNAME, as.character(DOBJ[[i]]), ", ")
  }
  DNAME = paste0(DNAME, "and ",as.character(DOBJ[[k]]))
  MNAME = riemobj$name
  # parameters
  myiter = max(50, round(maxiter))
  myeps  = min(max(as.double(eps),0),1e-5)  
  myperm = max(19, round(nperm))
  
  ## COMPUTATION
  #  Step 1. compute pooled frechet objective
  pooled.data = list()
  for (i in 1:k){
    pooled.data = c(pooled.data, datalist[[i]]$data)
  }
  pooled.ndata   = length(pooled.data)
  pooled.weight  = rep(1/pooled.ndata, pooled.ndata)
  frechet.pooled = inference_mean_intrinsic(MNAME, pooled.data, pooled.weight, myiter, myeps)
  
  # Step 2. compute individual frechet objective
  frechet.each = list()
  vec.ndata    = rep(0,k)
  for (i in 1:k){
    tgt.ndata    = length(datalist[[i]]$data)
    tgt.weight   = rep(1/tgt.ndata, tgt.ndata)
    vec.ndata[i] = tgt.ndata
    frechet.each[[i]] = inference_mean_intrinsic(MNAME, datalist[[i]]$data, tgt.weight, myiter, myeps)
  }
  
  # Step 3. get distance information only and compute via 'common_fanova'
  dist.pooled = as.vector(frechet.pooled$distvec)
  dist.class  = list()
  for (i in 1:k){
    dist.class[[i]] = as.vector(frechet.each[[i]]$distvec)
  }
  statinfo = common_fanova(dist.pooled, dist.class, MNAME, DNAME)
  Tnow     = as.double(statinfo$statistic)
  
  # Step 4. Monte Carlo simulation via permutation 
  Tvec = rep(0,myperm)
  for (i in 1:myperm){
    # permutation index generation
    randidx = aux_shuffle(vec.ndata)
    # compute individual-class statistic
    for (j in 1:k){
      tgt.data    = pooled.data[randidx[[j]]]
      tgt.ndata   = vec.ndata[j]
      tgt.weight  = rep(1/tgt.ndata, tgt.ndata)
      frechet.tmp = inference_mean_intrinsic(MNAME, tgt.data, tgt.weight, myiter, myeps)
      dist.class[[j]] = as.vector(frechet.tmp$distvec)
    }
    # compute the statistic
    statnow = common_fanova(dist.pooled, dist.class, MNAME, DNAME)
    Tvec[i] = as.double(statnow$statistic)
  }

  ############################################################
  # WRAP AND RETURN
  statinfo$p.value = (sum(Tvec >= Tnow)+1)/(myperm+1)
  return(statinfo)
}
#' @keywords internal
#' @noRd
common_fanova <- function(distall, distvecs, manifold, dataname){
  # get some parameters
  k = length(distvecs)  # number of classes
  n = length(distall)
  
  # get local information
  vec.nj    = rep(0,k)
  vec.Vj    = rep(0,k)
  vec.sig2j = rep(0,k)
  for (j in 1:k){
    distj = as.vector(distvecs[[j]])
    nj    = length(distj)
    
    vec.nj[j]    = nj
    vec.Vj[j]    = sum(distj^2)/nj
    vec.sig2j[j] = (sum(distj^4)/nj) - ((sum(distj^2)/nj)^2)
  }
  
  # get global information
  Vp   = sum(distall^2)/n
  lbdj = vec.nj/n
  
  # compute statistics
  Fn = Vp - sum(lbdj*vec.Vj)
  Un = 0
  for (j in 1:(k-1)){
    for (l in (j+1):k){
      Un = Un + ((lbdj[j]*lbdj[l])/(vec.sig2j[j]*vec.sig2j[l]))*((vec.Vj[j]-vec.Vj[l])^2)
    }
  }
  term1   = (n*Un)/sum(vec.nj/vec.sig2j)
  term2   = (n*(Fn^2))/sum((vec.nj^2)*vec.sig2j)
  thestat = term1+term2
  
  # compute p-value
  pvalue = stats::pchisq(thestat, df=(k-1), lower.tail = FALSE)
  
  # return output
  mfdname  = wrap_mfd2full(manifold)
  hname    = paste0("Frechet Analysis of Variance on ",mfdname," Manifold")
  Ha       = "at least one of equalities does not hold."
  names(thestat) = "Tn"
  
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name=dataname)
  class(res) = "htest"
  return(res)
}
