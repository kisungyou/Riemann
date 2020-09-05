## other functions for auxiliary
#  (1) wrap_vec2mat  : if a vector is given, return a matrix
#  (2) wrap_mfd2full : convert a simple manifold name into full expression
#  (3) aux_shuffle   : random permutation 
#  (4) aux_rmat2array3d & aux_rvec2array3d : riemdata object to 3d array

# (1) wrap_vec2mat : if a vector is given, return a matrix =====================
#' @keywords internal
#' @noRd
wrap_vec2mat <- function(x){
  if (is.vector(x)){
    return(matrix(x, ncol=1))
  } else if (is.matrix(x)){
    if ((nrow(x)==1)||(ncol(x)==1)){
      return(matrix(x, ncol=1))
    } else {
      return(x)
    }
  } else {
    stop("* wrap_vec2mat : internal error.")
  }
}
#  (2) wrap_mfd2full : convert a simple manifold name into full expression =====
#' @keywords internal
#' @noRd
wrap_mfd2full <- function(x){
  if (all(x=="euclidean")){
    return("Euclidean")
  } else if (all(x=="spd")){
    return("SPD")
  } else if (all(x=="correlation")){
    return("Correlation")
  } else if (all(x=="stiefel")){
    return("Stiefel")
  } else if (all(x=="grassmann")){
    return("Grassmann")
  } else if (all(x=="rotation")){
    return("Rotation")
  } else if (all(x=="multinomial")){
    return("Multinomial")
  } else if (all(x=="spdk")){
    return("SPD-k")
  } else if (all(x=="gaussian")){
    return("Gaussian")
  } else if (all(x=="sphere")){
    return("Sphere")
  } else {
    stop("* no such manifold name.")
  }
}

# (3) aux_shuffle : random permutation =========================================
#' @keywords internal
#' @noRd
aux_shuffle <- function(vecn){
  k = length(vecn)
  n = sum(vecn)
  allvec = seq_len(n)
  listvec = list()
  for (i in 1:(k-1)){
    vecnow = sample(allvec, vecn[i], replace = FALSE)
    listvec[[i]] = vecnow
    allvec = base::setdiff(allvec, vecnow)
  }
  listvec[[k]] = allvec
  return(listvec)
}

# (4) aux_rmat2array3d & aux_rvec2array3d : riemdata object to 3d array ========
#' @keywords internal
#' @noRd
aux_rmat2array3d <- function(riemobj){
  p = riemobj$size[1]
  r = riemobj$size[2]
  n = length(riemobj$data)
  
  output = array(0,c(p,r,n))
  for (i in 1:n){
    output[,,i] = riemobj$data[[i]]
  }
  return(output)
}
#' @keywords internal
#' @noRd
aux_rvec2array3d <- function(riemobj){
  p = length(riemobj$data[[1]])
  n = length(riemobj$data)
  
  output = array(0,c(p,1,n))
  for (i in 1:n){
    output[,,i] = matrix(riemobj$data[[i]], ncol=1)
  }
  return(output)
}