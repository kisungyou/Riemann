#' Prepare Data on SPD Manifold of Fixed-Rank
#' 
#' When \eqn{(p\times p)} SPD matrices are of fixed-rank \eqn{k < p}, they form 
#' a geometric structure represented by \eqn{(p\times k)} matrices,
#' \deqn{SPD(k,p) = \lbrace X \in \mathbf{R}^{(p\times p)}~\vert~ Y Y^\top = X, \textrm{rank}(X) = k \rbrace}
#' It's key difference from \eqn{\mathcal{S}_{++}^p} is that all matrices should be 
#' of fixed rank \eqn{k} where \eqn{k} is usually smaller than \eqn{p}. Inputs are 
#' given as \eqn{(p\times p)} matrices with specified \eqn{k} and \code{wrap.spdk} 
#' automatically decomposes input square matrices into rank-\eqn{k} representation matrices.
#' 
#' @param input data matrices to be wrapped as \code{riemdata} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a rank-\eqn{k} matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} matrices of rank-\eqn{k}.}
#' }
#' @param k rank of the SPD matrices.
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of \eqn{(p\times k)} representation of the corresponding rank-\eqn{k} SPSD matrices.}
#'   \item{size}{size of each representation matrix.}
#'   \item{name}{name of the manifold of interests, \emph{"spdk"}}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #                 Checker for Two Types of Inputs
#' #-------------------------------------------------------------------
#' #  Data Generation
#' d1 = array(0,c(10,10,3))
#' d2 = list()
#' for (i in 1:3){
#'   dat = matrix(rnorm(10*10),ncol=10)
#'   d1[,,i] = stats::cov(dat)
#'   d2[[i]] = d1[,,i]
#' }
#' 
#' #  Run
#' test1 = wrap.spdk(d1, k=2)
#' test2 = wrap.spdk(d2, k=2)
#' 
#' 
#' @references 
#' \insertRef{journee_lowrank_2010}{Riemann}
#' 
#' @concept wrapper
#' @export
wrap.spdk <- function(input, k){
  ## TAKE EITHER 3D ARRAY OR A LIST
  #  1. data format
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=FALSE)){
      stop("* wrap.spdk : input does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* wrap.spdk : input should be either a 3d array or a list.")
  }
  #  2. check all same size
  if (!check_list_eqsize(tmpdata, check.square=TRUE)){
    stop("* wrap.spdk : elements are not of same size.")
  }
  #  3. check and transform
  N = length(tmpdata)
  K = round(k)
  if ((K<1)||(K>nrow(tmpdata[[1]]))){
    stop("* wrap.spdk : target rank 'k' should be in [1,p]. For two extreme cases, use other geometries.")
  }
  for (n in 1:N){
    tmpdata[[n]] = single_spdkcheck(tmpdata[[n]], n, K)
  }
  
  ## WRAP AND RETURN THE S3 CLASS
  output = list()
  output$data = tmpdata
  output$size = dim(tmpdata[[1]])
  output$name = "spdk"
  return(structure(output, class="riemdata")) 
}
#' @keywords internal
#' @noRd
single_spdkcheck <- function(x, id, k){
  p = nrow(x)
  if (!(nrow(x)==ncol(x))){
    stop(paste0("* wrap.spdk : ",id,"-th element is not a square matrix."))
  }
  if (!isSymmetric(x)){
    stop(paste0("* wrap.spdk : ",id,"-th element is not a symmetric matrix."))
  }
  xrank = round(mat_rank(x))
  eigx  = base::eigen(x)
  if (xrank >= k){
    output = eigx$vectors[,1:k]%*%sqrt(diag(eigx$values[1:k]))
  } else {
    stop(paste0("* wrap.spdk : ",id,"-th element is rank deficient."))
  }
  return(output)
}
