# check_weight      : nonnegative numbers that sum to 1 of given length
# check_list_eqsize : for a list, all elements are of same size
# check_3darray     : check if 3d array of (p,p,N) type
# check_inputmfd    : check the object to abide by the structure
# check_spdmat      : check SPD matrix
# check_num_nonneg  : check a nonnegative real number
# check_unitvec     : check a unit-norm vector
# check_tworiems    : check whether two input 'riemdata' class are identical

# check_spdmat ------------------------------------------------------------
#' @keywords internal
#' @noRd
check_spdmat <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) < 1L || nrow(x) != ncol(x) || any(!is.finite(x))) return(FALSE)
  scale <- max(abs(x))
  if (scale == 0 || max(abs(x - t(x))) > 64 * .Machine$double.eps * scale) return(FALSE)
  !inherits(tryCatch(chol(x / 2 + t(x) / 2), error = identity), "error")
}

# check_weight      : nonnegative numbers that sum to 1 of given length ========
#' @keywords internal
#' @noRd
check_weight <- function(weight, N, fname) {
  if (!is.numeric(weight) || is.complex(weight) || !is.null(dim(weight)) ||
      length(weight) != N || any(!is.finite(weight)) ||
      any(weight < 0) || !any(weight > 0)) {
    stop("* ", fname, " : weights must be finite nonnegative numbers of length ",
         N, " with positive sum.", call. = FALSE)
  }
  scaled <- weight / max(weight)
  scaled / sum(scaled)
}

# check_3darray     : check if 3d array of (p,p,N) type ========================
#' @keywords internal
#' @noRd
check_3darray <- function(x, symmcheck = TRUE) {
  is.array(x) && length(dim(x)) == 3L && all(dim(x) > 0L) &&
    (!symmcheck || dim(x)[1L] == dim(x)[2L])
}

#' @keywords internal
#' @noRd
check_list_eqsize <- function(dlist, check.square = FALSE) {
  if (!is.list(dlist) || !length(dlist)) return(FALSE)
  first <- dlist[[1L]]
  if (is.null(dim(first))) {
    return(!check.square && length(first) > 0L &&
      all(vapply(dlist, function(x) is.null(dim(x)) && length(x) == length(first), logical(1))))
  }
  is.matrix(first) && all(dim(first) > 0L) &&
    (!check.square || nrow(first) == ncol(first)) &&
    all(vapply(dlist, function(x) is.matrix(x) && identical(dim(x), dim(first)), logical(1)))
}

# check_inputmfd    : check the object to abide by the structure ===============
#' @keywords internal
#' @noRd
check_inputmfd <- function(riemobj, funcname){
  mfdtype = strsplit(funcname,"[.]")[[1]][1]
  cond1 = (inherits(riemobj, "riemdata"))
  cond2 = all(riemobj$name==mfdtype)
  if (!(cond1&&cond2)){
    stop(paste0("* ",funcname," : input should be an object of 'riemdata' class with ",mfdtype,"-valued data."))
  }
}

# check_num_nonneg  : check a nonnegative real number ---------------------
#' @keywords internal
#' @noRd
check_num_nonneg <- function(x, funcname){
  cond1 = (length(x)==1)
  cond2 = ((all(is.finite(x)))&&(!any(is.na(x)))&&(all(x>=0)))
  if (cond1&&cond2){
    return(as.double(x))
  } else {
    stop(paste0("* ",funcname," : ",deparse(substitute(x))," is not a nonnegative number."))
  }
}

# check_unitvec : check a unit-norm vector --------------------------------
check_unitvec <- function(x, funcname){
  cond1 = is.vector(x)
  cond2 = (abs(sum(x^2)-1) < sqrt(.Machine$double.eps))
  if (cond1&&cond2){
    return(as.vector(x))
  } else {
    stop(paste0("* ",funcname," : ",deparse(substitute(x))," is not a unit vector."))
  }
}


# check_tworiems : check whether two input 'riemdata' class are id --------
#' @keywords internal
#' @noRd
check_tworiems <- function(riem1, riem2){
  # both are 'riemdata' classes
  if (!inherits(riem1,"riemdata")){
    return(FALSE)
  }
  if (!inherits(riem2,"riemdata")){
    return(FALSE)
  }
  # same manifold
  if (!all(riem1$name==riem2$name)){
    return(FALSE)
  }
  # same dimensionality
  if (!all(riem1$size==riem2$size)){
    return(FALSE)
  }
  return(TRUE)
}
