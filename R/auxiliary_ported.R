## auxiliary : ported functions from other packages 
#
#  T4transport_wassersteinD
#  T4transport_ipotD



# T4transport_ipotD -------------------------------------------------------
#' @keywords internal
#' @noRd
T4transport_ipotD <- function(D, p, wx, wy, lambda=1, ...){
  par_wx  = as.vector(wx)
  par_wy  = as.vector(wy)
  par_p   = max(1, as.double(p))
  par_lbd = max(sqrt(.Machine$double.eps), as.double(lambda))
  par_D   = as.matrix(D) 
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  par_iter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol  = max(sqrt(.Machine$double.eps), as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  par_inner = max(1, round(ifelse(("L"%in%pnames), params$L, 1)))

  output = cpp_ipot20(par_wx, par_wy, par_D, par_lbd, par_p, par_iter, par_tol, par_inner)
  return(output)
}
