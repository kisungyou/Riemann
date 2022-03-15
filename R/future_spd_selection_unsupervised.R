#' Unsupervised Variable Selection : Forward Type
#' 
#' 
#' 
#' 
#' 
#' @keywords internal
#' @noRd
spd.fusel <- function(spdobj, num.var=2){
  # PREPARE
  DNAME = paste0("'",deparse(substitute(spdobj)),"'")
  if ((!inherits(spdobj,"riemdata"))||(!all(spdobj$name=="spd"))){
    stop(paste0("* spd.fusel : input ",DNAME," should be an object of 'riemdata' class on 'spd' manifold.."))
  }
  num_var = max(1, round(num.var))
  p = round(spdobj$size[1])
  n = length(spdobj$data)
  if (num_var >= p){
    stop(paste0("* spd.fusel : 'num.var' should be smaller than ",p,"."))
  }
  
  # DATA WRANGLING + FRECHET MEAN COMPUTATION
  data3d = spd.wrap3d(spdobj$data) # (p x p x n)
  weight = rep(1/n, n)
  fmean  = inference_mean_extrinsic(spdobj$name, spdobj$data, weight, 100, 1e-6)$mean
  
  # STEP 1 : Select an entry of the largest variation
  id_all  = 1:p
  id_now  = spd_fusel_single(data3d, fmean)
  id_rest = setdiff(id_all, id_now)
  
  # STEP 2 : Iteratively
  if (num_var > 1){
    for (it in 2:num_var){
      tgt_length = length(id_rest)
      tgt_score  = rep(0, tgt_length)
      for (i in 1:tgt_length){
        tmp_id   = c(id_now, id_rest[i])
        tmp_data = data3d[tmp_id, tmp_id,]
        tmp_mean = fmean[tmp_id, tmp_id]
        tgt_score[i] = src_spd_variation(tmp_data, tmp_mean)
      }
      id_new  = id_rest[which.max(tgt_score)]
      id_now  = c(id_now, id_new)
      id_rest = setdiff(id_rest, id_new)
    }
  }
  
  # STEP 3 : Auxiliary Quantities
  var_full = src_spd_variation(data3d, fmean)/n
  if (num_var > 1){
    var_part = src_spd_variation(data3d[id_now,id_now,], fmean[id_now,id_now])/n
  } else {
    tgt_data = as.vector(data3d[id_now, id_now,])
    tgt_mean = as.double(fmean[id_now,id_now])
    var_part = sum((log(tgt_data) - log(tgt_mean))^2)/n
  }
  
  # RETURN
  output = list()
  output$index = id_now
  output$variation_full = var_full
  output$variation_part = var_part
  return(output)
}


# auxiliary function ------------------------------------------------------
#' @keywords internal
#' @noRd
spd_fusel_single <- function(data3d, fmean){
  p = base::nrow(fmean)
  score = rep(0, p)
  for (i in 1:p){
    tgt_data = as.vector(data3d[i,i,])
    tgt_mean = fmean[i,i]
    score[i] = sum((log(tgt_data) - log(tgt_mean))^2)
  }
  return(which.max(score))
}

# # personal example
# spd_mats = array(0,c(5,5,20))
# for (i in 1:10){
#   spd_mats[,,i] = stats::cov(matrix(rnorm(50*5), ncol=5))
# }
# for (j in 11:20){
#   randvec = stats::rnorm(50, sd=3)
#   randmat = cbind(sin(randvec), cos(randvec), sin(randvec)*cos(randvec))
#   randmat = cbind(randmat, matrix(rnorm(50*2), ncol=2))
#   spd_mats[,,j] = stats::cov(randmat)
# }
# spd_obj = wrap.spd(spd_mats)
# 
# spd.fusel(spd_obj, num.var = 1)
# spd.fusel(spd_obj, num.var = 2)
# spd.fusel(spd_obj, num.var = 3)
# spd.fusel(spd_obj, num.var = 4)
