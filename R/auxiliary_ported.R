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


# T4transport_wassersteinD ------------------------------------------------
#' @keywords internal
#' @noRd
T4transport_wassersteinD <- function(D, p, wx, wy){
  ## PARAMETERS 
  dxy = as.matrix(D)
  p   = max(1, double(p))
  wx  = as.vector(wx)
  wy  = as.vector(wy)
  
  cxy = (dxy^p)
  m   = length(wx); ww_m = matrix(wx, ncol=1)
  n   = length(wy); ww_n = matrix(wy, nrow=1)
  ones_m = matrix(rep(1,n),ncol=1)
  ones_n = matrix(rep(1,m),nrow=1)
  plan   = CVXR::Variable(m,n)
  
  wd.obj    <- CVXR::Minimize(CVXR::matrix_trace(t(cxy)%*%plan))
  wd.const1 <- list(plan >= 0)
  wd.const2 <- list(plan%*%ones_m==ww_m, ones_n%*%plan==ww_n)
  wd.prob   <- CVXR::Problem(wd.obj, c(wd.const1, wd.const2))
  wd.solve  <- CVXR::solve(wd.prob, solver="OSQP")
  
  if (all(wd.solve$status=="optimal")){ # successful
    gamma <- wd.solve$getValue(plan)
    value <- (base::sum(gamma*cxy)^(1/p))
    
    return(list(distance=value, plan=gamma))
  } else {                              # failed : use lpsolve
    cxy = (dxy^p)
    m   = nrow(cxy)
    n   = ncol(cxy)
    
    c  = as.vector(cxy)
    A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
    A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
    A  = rbind(A1, A2)
    
    f.obj = c
    f.con = A
    f.dir = rep("==",nrow(A))
    f.rhs = c(rep(1/m,m),rep(1/n,n))
    f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
    
    gamma = matrix(f.sol$solution, nrow=m)
    value = (sum(gamma*cxy)^(1/p))
    
    return(list(distance=value, plan=gamma))
  }
}
