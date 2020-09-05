# check_weight      : nonnegative numbers that sum to 1 of given length
# check_list_eqsize : for a list, all elements are of same size
# check_3darray     : check if 3d array of (p,p,N) type
# check_inputmfd    : check the object to abide by the structure


# check_weight      : nonnegative numbers that sum to 1 of given length ========
#' @keywords internal
#' @noRd
check_weight <- function(weight, N, fname){
  if ((!is.vector(weight))||(length(weight)!=N)){
    stop(paste0("* ",fname," : a weight parameter should be a vector of length corresponding to the provided data."))
  }
  if (any(weight)<= 0){
    stop(paste0("* ",fname," : we recommend to provide a weight vector of nonnegative weights."))
  }
  return(weight/base::sum(weight))
}

# check_3darray     : check if 3d array of (p,p,N) type ========================
#' @keywords internal
#' @noRd
check_3darray <- function(x, symmcheck=TRUE){
  cond1 = is.array(x)
  cond2 = (length(dim(x))==3)
  if (symmcheck){
    cond3 = (dim(x)[1] == dim(x)[2])  
  } else {
    cond3 = TRUE
  }
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# check_list_eqsize : for a list, all elements are of same size ================
#' @keywords internal
#' @noRd
check_list_eqsize <- function(dlist, check.square=FALSE){
  if (is.vector(dlist[[1]])){
    cond0 = all(unlist(lapply(dlist, is.vector))==TRUE)        # all vectors
    cond1 = (length(unique(unlist(lapply(dlist, length))))==1) # same length
    if (cond0&&cond1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    cond0 = all(unlist(lapply(dlist, is.matrix))==TRUE)      # all matrices
    cond1 = (length(unique(unlist(lapply(dlist, nrow))))==1) # same row size
    cond2 = (length(unique(unlist(lapply(dlist, ncol))))==1) # same col size
    if (check.square){
      cond3 = (nrow(dlist[[1]])==ncol(dlist[[1]]))
    } else {
      cond3 = TRUE
    }
    if (cond0&&cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
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