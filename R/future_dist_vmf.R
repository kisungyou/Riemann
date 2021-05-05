#' #' von Mises-Fisher Distribution
#' #' 
#' #' 
#' #' 
#' #' 
#' #' @name vmf
#' #' @concept distribution
#' #' @rdname vmf
#' NULL
#' 
#' #' @rdname vmf
#' #' @export
#' dvmf <- function(datalist, mu, kappa, log=FALSE){
#'   ## INITIALIZATION
#'   FNAME = "dvmf"
#'   myobj = wrap.sphere(datalist)
#'   mymu  = check_unitvec(mu, FNAME)
#'   mykap = check_num_nonneg(kappa, FNAME)
#'   
#' }
#' 
#' 
#' 
#' # auxiliary functions -----------------------------------------------------
#' #' @keywords internal
#' #' @noRd
#' rvmf_uniform <- function(n, mu, k=0){
#'   d = length(mu)
#'   x1 = matrix(rnorm(n*d),nrow=n)
#'   x  = x1/sqrt(base::rowSums(x1^2))
#'   return(x)
#' }
