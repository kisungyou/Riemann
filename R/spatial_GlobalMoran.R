#' Moran's \emph{I}
#' 
#' 
#' 
#' @keywords internal
#' @noRd
riem.moran <- function(riemobj, W, alternative=c("two.sided","greater","less"), ...){
  # ----------------------------------------------------------------------------
  # INPUTS : EXPLICIT
  # the data object
  if (!inherits(riemobj,"riemdata")){
    stop(paste0("* riem.moran : 'riemobj' should be an object of 'riemdata' class."))
  }
  N = length(riemobj$data)
  
  # weight
  if (!spatial_check_W(W, N)){
    stop("* riem.moran : 'W' is not a valid weight matrix. Please see the documentation.")
  }
  
  # hypothesis testing arguments
  par_H0    = match.arg(alternative)

  # INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("geometry"%in%pnames){
    par_geom = match.arg(tolower(params$geometry), c("intrinsic","extrinsic"))
  } else {
    par_geom = "intrinsic"
  }
  
  if ("ntest"%in%pnames){
    par_ntest = max(9, round(params$ntest))
  } else {
    par_ntest = 999
  }
  
  # ----------------------------------------------------------------------------
  # COMPUTE
  # weight normalize with zero diagonal
  diag(W) = 0 
  for (n in 1:N){
    tgt_W = as.vector(W[n,])
    W[n,] = tgt_W/base::sum(tgt_W)
  }
  
  # use CPP-accerelated routine
  cpprun = spatial_moran_global(riemobj$name,
                                par_geom,
                                riemobj$data,
                                W,
                                par_ntest)
  
  # ----------------------------------------------------------------------------
  # WRAP-UP
  out_statistic = as.double(cpprun$statistic)
  out_permuted  = as.vector(cpprun$permuted)
  
  if (identical(par_H0, "greater")){
    out_pvalue = (sum(out_statistic <= out_permuted)+1)/(length(out_permuted) + 1)
  } else if (identical(par_H0, "less")){
    out_pvalue = (sum(out_statistic >= out_permuted)+1)/(length(out_permuted) + 1)
  } else {
    pval1 = (sum(out_statistic <= out_permuted)+1)/(length(out_permuted) + 1)
    pval2 = (sum(out_statistic >= out_permuted)+1)/(length(out_permuted) + 1)
    out_pvalue = min(pval1,pval2)*2.0
  }
  
  # ----------------------------------------------------------------------------
  # RETURN
  output = list()
  output$statistic = as.double(cpprun$statistic)
  output$permuted  = as.vector(cpprun$permuted)
  output$p.value   = out_pvalue
  return(output)
}



# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
spatial_check_W <- function(W, n){
  cond1 = is.matrix(W)
  cond2 = all(W>=0)
  cond3 = ((base::nrow(W)==n)&&(base::ncol(W)==n))
  cond4 = all(is.finite(W))
  if (cond1&&cond2&&cond3&&cond4){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# ## personal test : comparison against the example from 'tmap' vignette.
# ## source - https://mgimond.github.io/simple_moransI_example/
# 
# library(sp)
# library(spdep)
# library(tmap)
# 
# # read the data
# s <- readRDS(url("https://github.com/mgimond/Data/raw/gh-pages/Exercises/fl_hr80.rds"))
# 
# # Define the neighbors (use queen case)
# nb <- poly2nb(s, queen=TRUE)
# 
# # Compute the neighboring average homicide rates
# lw <- nb2listw(nb, style="W", zero.policy=TRUE)
# 
# # Run the MC simulation version of the Moran's I test
# M1 <- moran.mc(s$HR80, lw, nsim=9999, alternative="less")
# 
# 
# # my way
# dat.euc = list()
# for (i in 1:length(s$HR80)){
#   dat.euc[[i]] = s$HR80[i]
# }
# dat.riem = wrap.euclidean(dat.euc)
# dat.W    = nb2mat(nb, zero.policy = TRUE, style="W")
# 
# # run my moran
# run_mine = riem.moran(dat.riem, dat.W, ntest=9999, alternative="less")
# 
# # visualize
# par(mfrow=c(1,2))
# plot(M1)
# 
# kdmine = density.default(run_mine$permuted)
# plot(kdmine)
# abline(v=run_mine$statistic, lwd=2, col="red")
