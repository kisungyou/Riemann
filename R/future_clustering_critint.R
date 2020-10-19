#' Internal Clustering Criteria
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @references 
#' \insertRef{dunn_fuzzy_1973a}{Riemann}
#' 
#' \insertRef{calinski_dendrite_1974}{Riemann}
#' 
#' \insertRef{baker_measuring_1975}{Riemann}
#' 
#' \insertRef{hubert_general_1976}{Riemann}
#' 
#' \insertRef{davies_cluster_1979}{Riemann}
#' 
#' \insertRef{rousseeuw_silhouettes_1987}{Riemann}
#' 
#' \insertRef{bezdek_new_1998}{Riemann}
#' 
#' \insertRef{saitta_bounded_2007}{Riemann}
#' 
#' 
#' @concept clustering
#' @keywords internal
#' @noRd
riem.critint <- function(riemobj, label, criteria, geometry=c("intrinsic","extrinsic")){
  ## PREPARE
  DNAME = paste0("'",deparse(substitute(riemobj)),"'") 
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.critint : input ",DNAME," should be an object of 'riemdata' class."))
  }
  nndata  = length(riemobj$data)
  mylabel = as.integer(as.factor(label))
  if (length(mylabel)!=nndata){
    stop(paste0("* riem.critint : input label should have length ",nndata,", which is the size of the given data."))
  }
  mygeom = ifelse(missing(geometry),"intrinsic",
                  match.arg(tolower(geometry),c("intrinsic","extrinsic")))
  
  ## RUN THE METHOD
  ncrit  = length(criteria)
  outval = list()
  for (i in 1:ncrit){
    mname = tolower(criteria[i])
    outval[[mname]] = switch(mname,
                             "dunn"  = cvi_internal_dunn(riemobj$name, mygeom, riemobj$data, mylabel),
                             "ch"    = cvi_internal_ch(riemobj$name, mygeom, riemobj$data, mylabel),
                             "ci"    = cvi_internal_ci(riemobj$name, mygeom, riemobj$data, mylabel),
                             "db"    = cvi_internal_db(riemobj$name, mygeom, riemobj$data, mylabel),
                             "sil"   = cvi_internal_sil(riemobj$name, mygeom, riemobj$data, mylabel),
                             "gamma" = cvi_internal_gamma(riemobj$name, mygeom, riemobj$data, mylabel),
                             "gd41"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 4, 1),
                             "gd43"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 4, 3),
                             "gd31"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 3, 1),
                             "gd33"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 3, 3),
                             "gd51"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 5, 1),
                             "gd53"  = cvi_internal_gdxx(riemobj$name, mygeom, riemobj$data, mylabel, 5, 3),
                             "score" = cvi_internal_score(riemobj$name, mygeom, riemobj$data, mylabel)
    )
  }
  return(outval)
}

# routines written in R ---------------------------------------------------
#' @keywords internal
#' @noRd
cvi_internal_sil <- function(mfd, geo, data, label){
  distobj     = stats::as.dist(basic_pdist(mfd, data, geo))
  func.import = utils::getFromNamespace("hidden_silhouette", "maotai")
  obj.silvals = func.import(distobj, label)
  return(obj.silvals$global)
}
#' @keywords internal
#' @noRd
cvi_internal_gamma <- function(mfd, geo, data, label){
  clindex = cvi_helper_R_classindex(label)
  distmat = basic_pdist(mfd, data, geo)
  diag(distmat) = Inf
  Nw = cvi_helper_R_nw(label)
  N  = length(label)
  K  = length(clindex)
  allid = 1:N
  
  val_den = 0
  for (k in 1:K){
    id_tgt  = clindex[[k]]
    id_else = setdiff(allid, id_tgt)
    ntgt = length(id_tgt)
    
    others = (distmat[id_else, id_else])
    others = as.vector(others[upper.tri(others)]) # only counting the pairs
    if (ntgt > 1){
      for (i in 1:(ntgt-1)){
        for (j in (i+1):ntgt){
          val_now = round(length(which(others < as.double(distmat[id_tgt[i], id_tgt[j]])))) 
          val_den = val_den + val_now
        }
      }
    }
  }
  
  val_num = ((N*(N-1)/2)-Nw)*Nw
  return(val_den/val_num)
}


# extra routines for writing in R -----------------------------------------
#' @keywords internal
#' @noRd
cvi_helper_R_nw <- function(label){
  clindex = cvi_helper_R_classindex(label)
  K = length(clindex)
  
  output = 0
  for (k in 1:K){
    ck = length(clindex[[k]])
    output = output + (ck*(ck-1)/2)
  }
  return(output)
}
#' @keywords internal
#' @noRd
cvi_helper_R_classindex <- function(label){
  K = length(unique(label))
  output = list()
  for (k in 1:K){
    output[[k]] = which(label==k)
  }
  return(output)
}
#' @keywords internal
#' @noRd
cvi_helper_R_classmean <- function(mfd, geo, data, label){
  classindex = cvi_helper_R_classindex(label)
  
  output = list()
  K      = length(classindex)
  for (k in 1:K){
    partdat   = data[classindex[[k]]]
    partnum   = length(partdat)
    tmpweight = rep(1/partnum, partnum)
    if (partnum < 2){
      output[[k]] = partdat[[1]]
    } else {
      if (all(geo=="intrinsic")){
        runk = inference_mean_intrinsic(mfd, partdat, tmpweight, 100, 1e-5)  
      } else {
        runk = inference_mean_extrinsic(mfd, partdat, tmpweight, 100, 1e-5)  
      }
      output[[k]] = runk$mean
    }
  }
  return(output)
}





# ## GENERATE DATA
# mydata = list()
# for (i in 1:10){
#   tgt = c(1, stats::rnorm(2, sd=0.1))
#   mydata[[i]] = tgt/sqrt(sum(tgt^2))
# }
# for (i in 11:20){
#   tgt = c(rnorm(1,sd=0.1),1,rnorm(1,sd=0.1))
#   mydata[[i]] = tgt/sqrt(sum(tgt^2))
# }
# for (i in 21:30){
#   tgt = c(stats::rnorm(2, sd=0.1), 1)
#   mydata[[i]] = tgt/sqrt(sum(tgt^2))
# }
# myriem = wrap.sphere(mydata)
# # 
# ## K-MEANS++ WITH K=2,3,4
# vec_x = 2:11
# vec_y = rep(0,10)
# for (i in 1:10){
#   cli = riem.kmedoids(myriem, k=(i+1))$cluster
#   vec_y[i] = riem.critint(myriem, cli, criteria = "gd31")$'gd31'
# }
# plot(vec_x, vec_y, "b")